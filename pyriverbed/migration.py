"""Meander migration and cutoffs.

The migration model is the closed-form solution of the linearised bend theory
of Ikeda, Parker & Sawai (1981): the dimensionless near-bank excess velocity
is the sum of a decaying inlet transient, a local-curvature term, and a term
proportional to the phase-lagged curvature. Bank retreat is taken proportional
to that velocity.

Two cutoff mechanisms are implemented, and they are different in kind:

* a **neck cutoff** happens when a loop grows until its two limbs touch, so it
  is detected geometrically and is deterministic;
* a **chute cutoff** happens when flow carves a new channel across the
  floodplain, so its initiation depends on events the model does not resolve
  and is treated as conditionally random.

See ``THEORY_GUIDE.md`` sections 9 and 10.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from .config import (ChannelConfig, ChuteCutoffConfig, MigrationConfig,
                     NeckCutoffConfig)
from .logging_utils import get_logger
from .planform import Centerline, arc_length, smooth

log = get_logger(__name__)

__all__ = ["Cutoff", "MigrationResult", "migrate", "find_neck_cutoff",
           "find_chute_entrances", "chute_valley_angle", "find_chute_cutoff",
           "carve_cutoff", "carve_channel", "end_taper"]


@dataclass
class Cutoff:
    """A cutoff event and the oxbow lake it left behind.

    Attributes
    ----------
    kind
        ``'neck'`` or ``'chute'``.
    step
        Time step at which it happened.
    entrance, exit
        Node indices bounding the bypassed reach, on the pre-cutoff centerline.
    oxbow_x, oxbow_y
        Coordinates of the abandoned reach.
    """

    kind: str
    step: int
    entrance: int
    exit: int
    oxbow_x: np.ndarray = field(default_factory=lambda: np.zeros(0))
    oxbow_y: np.ndarray = field(default_factory=lambda: np.zeros(0))

    @property
    def bypassed_nodes(self) -> int:
        """Number of nodes removed from the centerline."""
        return int(self.oxbow_x.size)


@dataclass
class MigrationResult:
    """Outcome of one migration time step.

    Attributes
    ----------
    centerline
        The migrated centerline.
    mean_rate
        Mean nodal displacement over the step (m).
    near_bank_velocity
        The dimensionless near-bank excess velocity field.
    """

    centerline: Centerline
    mean_rate: float
    near_bank_velocity: np.ndarray


def near_bank_velocity(line: Centerline, cur_filtered: np.ndarray,
                       cur_lagged: np.ndarray, channel: ChannelConfig,
                       migration: MigrationConfig,
                       rng: np.random.Generator) -> np.ndarray:
    """Dimensionless near-bank excess velocity, after Ikeda et al. (1981).

    .. math:: u_b(\\tilde{s}) = a_1 e^{-a_2 \\tilde{s}} + a_3 \\tilde{C}
              + a_4 \\tilde{C}^{lag}

    The local-curvature coefficient ``a_3`` is negative while the lagged
    coefficient ``a_4`` is positive and several times larger for realistic
    parameters. It is that phase-shifted term which makes meanders both grow
    and translate downstream; a model without the lag would not meander.

    ``chi``, the cube root of the inverse sinuosity, is the feedback of the
    planform on the hydraulics: as the river lengthens its slope falls, so its
    velocity scale falls with the cube root of the slope. This is what makes
    long runs self-limiting.

    Parameters
    ----------
    line
        The current centerline.
    cur_filtered, cur_lagged
        Filtered and phase-lagged curvature (1/m).
    channel, migration
        Configuration.
    rng
        Random generator, used for the inlet perturbation.

    Returns
    -------
    ndarray
        The near-bank excess velocity at every node.
    """
    width = channel.width
    beta = channel.beta
    scour = channel.scour_factor

    s_hat = line.s / width
    cur_hat = cur_filtered * width
    lag_hat = cur_lagged * width

    valley = line.valley_length
    chi = (valley / line.length) ** (1 / 3) if line.length > 0 else 1.0

    a1 = migration.ub0 * (2 * rng.random() - 1) + chi * migration.c0
    a2 = 2 * migration.cf0 * beta * chi
    a3 = -chi
    a4 = migration.cf0 * beta * (
        chi ** 5 * migration.fr0 ** 2
        + (scour + 1) * chi ** 2
        + 5 * np.sqrt(migration.cf0) * (scour + chi ** 2 * migration.fr0 ** 2)
    )
    return a1 * np.exp(-a2 * s_hat) + a3 * cur_hat + a4 * lag_hat


def end_taper(line: Centerline, channel: ChannelConfig,
              migration: MigrationConfig) -> np.ndarray:
    """Weights that ramp the bank displacement to zero at both reach ends.

    The linear theory has no upstream reach over which to spread an inlet
    perturbation, so at a free inlet the ``ub0`` noise and the curvature
    feedback amplify each other into a hook of curvature far beyond the range
    where the theory holds, which then coils up into spurious neck cutoffs.
    Fixing the end nodes is the standard remedy and leaves the interior
    untouched.

    The ramp is a raised cosine, so the weights and their first derivative are
    both continuous where the taper meets the interior; a linear ramp would put
    a weak kink there.

    Parameters
    ----------
    line
        The current centerline.
    channel
        Channel configuration, for the width scale.
    migration
        Migration configuration, for ``end_taper_widths``.

    Returns
    -------
    ndarray
        One weight in [0, 1] per node. All ones if the taper is switched off.
    """
    length = migration.end_taper_widths * channel.width
    if length <= 0 or line.length <= 0:
        return np.ones(line.n_nodes)
    # Never taper more than a third of the reach from each end.
    length = min(length, line.length / 3)
    distance = np.minimum(line.s, line.s[-1] - line.s)
    ramp = np.clip(distance / length, 0.0, 1.0)
    return 0.5 * (1.0 - np.cos(np.pi * ramp))


def migrate(line: Centerline, cur_filtered: np.ndarray,
            cur_lagged: np.ndarray, channel: ChannelConfig,
            migration: MigrationConfig,
            rng: np.random.Generator) -> MigrationResult:
    """Advance the centerline by one migration time step.

    Every node is displaced along the local normal by ``E0 * u_b * dt``, the
    Ikeda et al. (1981) bank erosion law. Both banks move together, so the
    channel width never changes.

    Parameters
    ----------
    line
        The current centerline.
    cur_filtered, cur_lagged
        Filtered and phase-lagged curvature (1/m).
    channel, migration
        Configuration.
    rng
        Random generator.

    Returns
    -------
    MigrationResult
    """
    width = channel.width
    ub = near_bank_velocity(line, cur_filtered, cur_lagged, channel,
                            migration, rng)

    displacement = migration.e0 * ub * migration.dt      # in channel widths
    displacement = displacement * end_taper(line, channel, migration)
    x = line.x / width + displacement * np.sin(line.theta)
    y = line.y / width - displacement * np.cos(line.theta)
    x, y = x * width, y * width

    rate = float(np.mean(np.hypot(x - line.x, y - line.y)))
    moved = Centerline(x=x, y=y, s=arc_length(x, y),
                       curvature=line.curvature.copy(), theta=line.theta.copy())
    return MigrationResult(centerline=moved, mean_rate=rate,
                           near_bank_velocity=ub)


# --------------------------------------------------------------------------
# neck cutoff
# --------------------------------------------------------------------------

def find_neck_cutoff(line: Centerline, channel: ChannelConfig,
                     config: NeckCutoffConfig) -> tuple[int, int]:
    """Find a neck cutoff: two nodes far apart along the channel but close in space.

    v1.x tested every node pair, an O(n^2) scan that dominated the cost of a
    long run. Here a k-d tree returns only the pairs that are actually within
    the threshold distance, which is O(n log n) and gives the same answer: the
    first qualifying pair in ascending node order.

    ``config.end_margin_widths`` of each end is excluded from the search. Those
    stretches carry the artificial straight extensions and the inlet transient
    of the migration model, so a pair found there reflects the boundary
    treatment rather than a meander loop closing on itself.

    Parameters
    ----------
    line
        The current centerline.
    channel
        Channel configuration, for the width scale.
    config
        Neck cutoff settings.

    Returns
    -------
    entrance, exit : int
        Node indices bounding the reach to bypass, or ``(-1, -1)`` if there is
        no neck cutoff.
    """
    if not config.enabled or line.n_nodes < 4:
        return -1, -1
    spacing = line.mean_spacing
    if spacing <= 0:
        return -1, -1
    min_gap = config.min_separation_widths * channel.width / spacing
    max_dist = config.max_distance_widths * channel.width

    points = np.column_stack((line.x, line.y))
    tree = cKDTree(points)
    pairs = tree.query_pairs(max_dist, output_type="ndarray")
    if pairs.size == 0:
        return -1, -1
    margin = config.end_margin_widths * channel.width
    first = int(np.searchsorted(line.s, margin))
    last = int(np.searchsorted(line.s, line.s[-1] - margin))

    lo = np.minimum(pairs[:, 0], pairs[:, 1])
    hi = np.maximum(pairs[:, 0], pairs[:, 1])
    keep = (hi - lo > min_gap) & (lo >= max(first, 1)) & (hi <= last)
    if not np.any(keep):
        return -1, -1
    lo, hi = lo[keep], hi[keep]
    # Match v1.x ordering: smallest entrance first, then smallest exit.
    order = np.lexsort((hi, lo))
    return int(lo[order[0]]), int(hi[order[0]])


# --------------------------------------------------------------------------
# chute cutoff
# --------------------------------------------------------------------------

def find_chute_entrances(cur_filtered: np.ndarray,
                         config: ChuteCutoffConfig) -> np.ndarray:
    """Locate the candidate entrance/exit points of chute channels.

    Inflection points are where the curvature changes sign. Bend apexes are
    the points of locally maximum absolute curvature between two consecutive
    inflection points. Which family is used is set by ``config.entrance``:
    ``'apex'`` is the bar-surface route, ``'inflection'`` cuts a whole bend out
    of the planform.

    Parameters
    ----------
    cur_filtered
        Filtered curvature.
    config
        Chute cutoff settings.

    Returns
    -------
    ndarray
        Node indices of the candidate entrance/exit points, ascending.
    """
    inflections = np.where(np.diff(np.sign(cur_filtered)) != 0)[0]
    if config.entrance == "inflection":
        return inflections
    if inflections.size < 2:
        return np.zeros(0, dtype=int)
    magnitude = np.abs(cur_filtered)
    apexes = np.empty(inflections.size - 1, dtype=int)
    for k in range(inflections.size - 1):
        start, end = inflections[k], inflections[k + 1]
        apexes[k] = start + int(np.argmax(magnitude[start:end]))
    return apexes


def chute_valley_angle(x: np.ndarray, y: np.ndarray, i: int, j: int) -> float:
    """Angle between a chute chord and the valley axis, in degrees.

    The chute channel is approximated by the straight chord from node *i* to
    node *j*, and the valley axis by the chord between the two ends of the
    centerline. A small angle means the chute runs down the valley and so
    enjoys the largest slope advantage over the along-channel path; a chute
    perpendicular to the valley has none.

    Measuring against the centerline's own end-to-end chord rather than
    against the x axis makes the test invariant to how the reach happens to be
    rotated in the coordinate system.

    Returns
    -------
    float
        The angle in ``[0, 90]`` degrees.
    """
    cx, cy = x[j] - x[i], y[j] - y[i]
    vx, vy = x[-1] - x[0], y[-1] - y[0]
    return float(np.degrees(np.arctan2(abs(cx * vy - cy * vx),
                                       abs(cx * vx + cy * vy))))


def find_chute_cutoff(line: Centerline, cur_filtered: np.ndarray,
                      step: int, channel: ChannelConfig,
                      config: ChuteCutoffConfig,
                      rng: np.random.Generator) -> tuple[int, int]:
    """Find a chute cutoff for the current time step.

    A chute cutoff is triggered with probability ``config.frequency`` per time
    step, after a spin-up of ``config.start_step`` steps. Once triggered, all
    admissible chute channels are collected and one is drawn uniformly. A
    candidate spanning ``config.span`` entrance points is admissible when the
    reach it bypasses is at least ``config.min_length_widths`` channel widths
    long, that reach is at least ``config.min_sinuosity`` times longer than the
    chute chord, the chord's angle with the valley axis is at most
    ``config.max_valley_angle`` degrees, and both ends stay
    ``config.end_margin`` entrance points clear of the reach ends.

    Parameters
    ----------
    line
        The current centerline.
    cur_filtered
        Filtered curvature, indexed consistently with *line*.
    step
        Current time step.
    channel, config
        Configuration.
    rng
        Random generator.

    Returns
    -------
    entrance, exit : int
        Node indices bounding the reach to bypass, or ``(-1, -1)`` if there is
        no chute cutoff this step.
    """
    if not config.enabled or step < config.start_step:
        return -1, -1
    if config.frequency <= 0 or rng.random() >= config.frequency:
        return -1, -1

    entrances = find_chute_entrances(cur_filtered, config)
    min_nodes = config.min_length_widths * 2 * channel.n_offsets
    first = config.end_margin
    last = entrances.size - config.span - config.end_margin
    candidates: list[tuple[int, int]] = []
    for k in range(first, last):
        i, j = int(entrances[k]), int(entrances[k + config.span])
        if j - i <= min_nodes:
            continue
        if chute_valley_angle(line.x, line.y, i, j) > config.max_valley_angle:
            continue
        # The slope advantage: how much shorter the chute is than the reach it
        # replaces. Without this test a chute can be carved across an almost
        # straight reach, which gains nothing and does not happen in nature.
        chord = float(np.hypot(line.x[j] - line.x[i], line.y[j] - line.y[i]))
        if chord <= 0:
            continue
        if (line.s[j] - line.s[i]) / chord < config.min_sinuosity:
            continue
        candidates.append((i, j))
    if not candidates:
        return -1, -1
    return candidates[int(rng.integers(len(candidates)))]


# --------------------------------------------------------------------------
# carving
# --------------------------------------------------------------------------

def carve_channel(x: np.ndarray, y: np.ndarray, index: int, spacing: float,
                  half_window: int, n_passes: int = 20
                  ) -> tuple[np.ndarray, np.ndarray]:
    """Turn the bare join left by a cutoff into an actual channel.

    Concatenating ``x[:i+1]`` with ``x[j:]`` leaves the whole cutoff channel as
    a *single* segment -- for a chute across a couple of meander loops, one
    segment several kilometres long between two nodes ten metres apart -- and a
    hard corner at each of its ends. Two things then go wrong: the corner is
    tighter than any channel can bend, and once the next resampling populates
    the long segment the corner shows up as a hairpin that the neck cutoff
    detector reads as a loop closing on itself. Every chute cutoff was
    spawning a spurious neck cutoff at its own junction on the following step.

    So the new channel is given nodes at the prevailing spacing, and then each
    of its two junctions is rounded by smoothing a window around it. Only
    then is the window uniformly sampled, which is what a Savitzky-Golay filter
    needs: applied to the bare join it would drag the neighbouring nodes
    hundreds of metres along the chord instead of rounding anything.

    Physically this is the young cutoff channel establishing itself and
    reworking its own junctions, which is where it should happen.

    Parameters
    ----------
    x, y
        Coordinates of the joined centerline, with the cutoff channel as the
        single segment from ``index`` to ``index + 1``.
    index
        Node index of the upstream end of the cutoff channel.
    spacing
        Target node spacing (m), normally the reach's mean spacing.
    half_window
        Nodes to smooth on each side of each junction.
    n_passes
        Savitzky-Golay passes to apply inside each window.

    Returns
    -------
    x, y : ndarray
        Coordinates with the cutoff channel discretised and its junctions
        rounded.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if index + 1 >= x.size or spacing <= 0:
        return x.copy(), y.copy()

    x0, y0, x1, y1 = x[index], y[index], x[index + 1], y[index + 1]
    chord = float(np.hypot(x1 - x0, y1 - y0))
    n_new = int(chord // spacing)
    if n_new > 1:
        fraction = np.arange(1, n_new)[:, None] / n_new
        inserted = np.array([x0, y0]) + fraction * np.array([x1 - x0, y1 - y0])
        x = np.concatenate((x[:index + 1], inserted[:, 0], x[index + 1:]))
        y = np.concatenate((y[:index + 1], inserted[:, 1], y[index + 1:]))
    else:
        n_new = 1

    x, y = x.copy(), y.copy()
    for junction in (index, index + n_new):
        lo = max(junction - half_window, 0)
        hi = min(junction + half_window + 1, x.size)
        if hi - lo < 5:                    # too short for a 5-point filter
            continue
        _, xs, ys = smooth(x[lo:hi], y[lo:hi], n_passes)
        x[lo:hi], y[lo:hi] = xs, ys
    return x, y


def carve_cutoff(line: Centerline, i: int, j: int, kind: str, step: int,
                 channel: ChannelConfig | None = None
                 ) -> tuple[Centerline, Cutoff]:
    """Remove the reach between nodes *i* and *j*, returning it as an oxbow lake.

    The new channel is discretised and its junctions rounded by
    :func:`carve_channel` over about one channel width to each side. Without
    that, the join is a single long segment ending in corners that the light
    per-step smoothing cannot remove, and they immediately trigger a further,
    spurious neck cutoff at the same place.

    Parameters
    ----------
    line
        The centerline to cut.
    i, j
        Node indices bounding the reach to bypass.
    kind
        ``'neck'`` or ``'chute'``, recorded on the returned :class:`Cutoff`.
    step
        Current time step.
    channel
        Channel configuration, used to size the rounding window in channel
        widths. ``None`` skips the rounding.

    Returns
    -------
    centerline : Centerline
        The shortened centerline.
    cutoff : Cutoff
        The event, carrying the oxbow lake.
    """
    oxbow_x = line.x[i + 1:j].copy()
    oxbow_y = line.y[i + 1:j].copy()
    x = np.concatenate((line.x[:i + 1], line.x[j:]))
    y = np.concatenate((line.y[:i + 1], line.y[j:]))
    if channel is not None:
        spacing = line.mean_spacing
        if spacing > 0:
            half = max(int(round(channel.width / spacing)), 2)
            x, y = carve_channel(x, y, i, spacing, half)
    cut = Centerline(x=x, y=y, s=arc_length(x, y),
                     curvature=np.zeros(x.size), theta=np.zeros(x.size))
    return cut, Cutoff(kind=kind, step=step, entrance=i, exit=j,
                       oxbow_x=oxbow_x, oxbow_y=oxbow_y)
