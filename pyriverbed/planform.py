"""Channel planform: generation, resampling, smoothing, curvature, phase lag.

Everything in this module works on a *centerline*, represented by the
:class:`Centerline` value object. Functions take and return plain arrays or
``Centerline`` instances and never touch global state, which is what lets the
same process run several configurations one after another -- something v1.x
could not do, because its Numba kernels froze the module-level parameters at
compile time.

All the numerics are vectorised NumPy. v1.x needed Numba to make its explicit
Python loops bearable; expressing the same operations as array operations is
both faster and one dependency lighter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter

from .config import ZERO, ChannelConfig, Config, CurvatureConfig, KinoshitaConfig
from .logging_utils import get_logger

log = get_logger(__name__)

__all__ = [
    "Centerline",
    "arc_length",
    "build_kinoshita",
    "curvature",
    "despike_curvature",
    "extend_ends",
    "filter_curvature",
    "kinoshita_equation_text",
    "load_centerline",
    "phase_lag",
    "resample",
    "smooth",
]


# --------------------------------------------------------------------------
# value object
# --------------------------------------------------------------------------

@dataclass
class Centerline:
    """A river centerline and the quantities derived from it.

    Attributes
    ----------
    x, y
        Node coordinates (m).
    s
        Cumulative arc length from the upstream end (m).
    curvature
        Signed curvature (1/m), positive for a left turn.
    theta
        Direction angle of the local tangent (rad).
    """

    x: np.ndarray
    y: np.ndarray
    s: np.ndarray
    curvature: np.ndarray
    theta: np.ndarray

    @property
    def n_nodes(self) -> int:
        """Number of nodes."""
        return int(self.x.size)

    @property
    def length(self) -> float:
        """Total centerline length (m)."""
        return float(self.s[-1])

    @property
    def valley_length(self) -> float:
        """Straight-line distance between the two ends (m)."""
        return float(np.hypot(self.x[-1] - self.x[0], self.y[-1] - self.y[0]))

    @property
    def sinuosity(self) -> float:
        """Ratio of channel length to valley length."""
        valley = self.valley_length
        return float(self.length / valley) if valley > 0 else float("inf")

    @property
    def mean_spacing(self) -> float:
        """Mean node spacing (m)."""
        return float(np.mean(np.diff(self.s))) if self.n_nodes > 1 else 0.0

    def copy(self) -> "Centerline":
        """Return a deep copy."""
        return Centerline(self.x.copy(), self.y.copy(), self.s.copy(),
                          self.curvature.copy(), self.theta.copy())

    def spacing_stats(self) -> tuple[float, float, float]:
        """Return the mean, median and mode of the node spacing (m).

        The mode is taken over spacings rounded to the centimetre, which is
        what makes it meaningful for hand-digitised centerlines.
        """
        d = np.diff(self.s)
        if d.size == 0:
            return 0.0, 0.0, 0.0
        rounded = np.round(d, 2)
        values, counts = np.unique(rounded, return_counts=True)
        return (float(np.mean(d)), float(np.median(d)),
                float(values[int(np.argmax(counts))]))


# --------------------------------------------------------------------------
# basic geometry
# --------------------------------------------------------------------------

def arc_length(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Cumulative arc length along the polyline ``(x, y)``, starting at 0."""
    steps = np.hypot(np.diff(x), np.diff(y))
    return np.concatenate(([0.0], np.cumsum(steps)))


def resample(x: np.ndarray, y: np.ndarray,
             spacing: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resample a polyline onto equally spaced nodes.

    The number of intervals is chosen so that the actual spacing is as close
    to *spacing* as possible without changing the total length, and both end
    nodes are reproduced exactly.

    Parameters
    ----------
    x, y
        Input node coordinates.
    spacing
        Target node spacing (m).

    Returns
    -------
    s, x, y : ndarray
        Arc length and coordinates of the resampled polyline.
    """
    s = arc_length(x, y)
    total = float(s[-1])
    if total <= 0:
        return s, np.array(x, dtype=float), np.array(y, dtype=float)
    n_intervals = int(total // spacing) + 1
    stations = np.linspace(0.0, total, n_intervals + 1)
    xn = np.interp(stations, s, x)
    yn = np.interp(stations, s, y)
    xn[0], yn[0], xn[-1], yn[-1] = x[0], y[0], x[-1], y[-1]
    return arc_length(xn, yn), xn, yn


def smooth(x: np.ndarray, y: np.ndarray,
           n_passes: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Smooth a polyline with repeated Savitzky-Golay passes.

    A 5-point, second-order Savitzky-Golay filter is applied to ``x`` and
    ``y`` independently, *n_passes* times. Unlike a moving average it
    preserves the position and amplitude of the bend apexes, which is what
    keeps the pool depths right. The two end nodes are restored after every
    pass so that smoothing cannot retract the reach.

    Parameters
    ----------
    x, y
        Input node coordinates.
    n_passes
        Number of passes. Zero returns the input unchanged.

    Returns
    -------
    s, x, y : ndarray
        Arc length and coordinates of the smoothed polyline.
    """
    x = np.asarray(x, dtype=float).copy()
    y = np.asarray(y, dtype=float).copy()
    if n_passes <= 0 or x.size < 5:
        return arc_length(x, y), x, y
    xa, xb, ya, yb = x[0], x[-1], y[0], y[-1]
    for _ in range(int(n_passes)):
        x = savgol_filter(x, 5, 2, mode="nearest")
        y = savgol_filter(y, 5, 2, mode="nearest")
        x[0], x[-1], y[0], y[-1] = xa, xb, ya, yb
    return arc_length(x, y), x, y


# --------------------------------------------------------------------------
# curvature
# --------------------------------------------------------------------------

def _segment_angles(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Direction angle of every segment of the polyline (rad)."""
    return np.arctan2(np.diff(y), np.diff(x))


def _tangent_angles(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Per-node tangent angle (rad), repeating the last segment at the end."""
    angles = _segment_angles(x, y)
    return np.concatenate((angles, angles[-1:]))


def curvature_arctan2(s: np.ndarray, x: np.ndarray,
                      y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Signed curvature from the change in tangent direction.

    The direction change between the segment leaving a node and the segment
    arriving at it, divided by the arc length it spans. This is pyRiverBed's
    default estimator.

    The angle difference is wrapped into ``(-pi, pi]``, which v1.x did not
    do. Without wrapping, a reach whose flow direction crosses due west picks
    up a spurious 2*pi jump and therefore a huge false curvature spike. The
    wrapping is a no-op for every well-resolved centerline.

    Returns
    -------
    curvature, theta : ndarray
        Signed curvature (1/m) and per-node tangent angle (rad).
    """
    n = x.size
    cur = np.zeros(n)
    if n < 3:
        return cur, _tangent_angles(x, y)
    angles = _segment_angles(x, y)
    delta = np.diff(angles)
    delta = (delta + np.pi) % (2 * np.pi) - np.pi
    span = s[2:] - s[:-2]
    with np.errstate(divide="ignore", invalid="ignore"):
        cur[1:-1] = np.where(span > 0, 2.0 * delta / span, 0.0)
    cur[np.abs(cur) < ZERO] = 0.0
    return cur, np.concatenate((angles, angles[-1:]))


def curvature_cosine(s: np.ndarray, x: np.ndarray,
                     y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Signed curvature via the law of cosines.

    Numerically a different route to the same quantity as
    :func:`curvature_arctan2`: the turning angle comes from ``arccos`` of the
    normalised dot product of consecutive segments, and its sign from their
    cross product. Provided for comparison.
    """
    n = x.size
    cur = np.zeros(n)
    if n < 3:
        return cur, _tangent_angles(x, y)
    dx, dy = np.diff(x), np.diff(y)
    norm = np.hypot(dx, dy)
    a_x, a_y = dx[1:], dy[1:]
    b_x, b_y = dx[:-1], dy[:-1]
    denom = norm[1:] * norm[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = np.where(denom > 0, (a_x * b_x + a_y * b_y) / denom, 1.0)
    angle = np.arccos(np.clip(cosine, -1.0, 1.0))
    angle = np.copysign(angle, a_x * b_y * -1.0 + a_y * b_x * 0.0
                        + (b_x * a_y - b_y * a_x))
    span = s[2:] - s[:-2]
    with np.errstate(divide="ignore", invalid="ignore"):
        cur[1:-1] = np.where(span > 0, 2.0 * angle / span, 0.0)
    cur[np.abs(cur) < ZERO] = 0.0
    return cur, _tangent_angles(x, y)


def curvature_circumcircle(s: np.ndarray, x: np.ndarray,
                           y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Signed curvature from the circle through three consecutive nodes.

    The curvature is ``1/R`` with ``R = abc / (4 * area)`` for the triangle
    formed by ``P[i-1], P[i], P[i+1]``. v1.x offered this estimator unsigned,
    which made it unusable for anything but a magnitude check; here the sign
    is taken from the cross product so it can actually be compared with the
    other two.
    """
    n = x.size
    cur = np.zeros(n)
    if n < 3:
        return cur, _tangent_angles(x, y)
    a = np.hypot(x[2:] - x[1:-1], y[2:] - y[1:-1])
    b = np.hypot(x[2:] - x[:-2], y[2:] - y[:-2])
    c = np.hypot(x[1:-1] - x[:-2], y[1:-1] - y[:-2])
    cross = ((x[1:-1] - x[:-2]) * (y[2:] - y[1:-1])
             - (y[1:-1] - y[:-2]) * (x[2:] - x[1:-1]))
    area = 0.5 * np.abs(cross)
    with np.errstate(divide="ignore", invalid="ignore"):
        radius = np.where(area > 0, a * b * c / (4.0 * area), np.inf)
        kappa = np.where(np.isfinite(radius) & (radius > 0), 1.0 / radius, 0.0)
    cur[1:-1] = np.copysign(kappa, cross)
    cur[np.abs(cur) < ZERO] = 0.0
    return cur, _tangent_angles(x, y)


_CURVATURE_METHODS = {
    "arctan2": curvature_arctan2,
    "cosine": curvature_cosine,
    "circumcircle": curvature_circumcircle,
}


def despike_curvature(cur: np.ndarray) -> np.ndarray:
    """Replace isolated curvature spikes by the mean of their neighbours.

    A node is a spike when it departs from the neighbour mean by more than
    five times the neighbour-to-neighbour difference. The sweep is sequential
    -- a corrected value is used when testing the next node -- so that a run
    of adjacent spikes is walked down rather than left half-fixed.
    """
    cur = np.asarray(cur, dtype=float).copy()
    for i in range(1, cur.size - 1):
        mean = 0.5 * (cur[i - 1] + cur[i + 1])
        if abs(cur[i] - mean) > 5.0 * abs(cur[i - 1] - cur[i + 1]):
            cur[i] = mean
    return cur


def curvature(s: np.ndarray, x: np.ndarray, y: np.ndarray,
              config: CurvatureConfig | None = None
              ) -> tuple[np.ndarray, np.ndarray]:
    """Estimate curvature using the configured method, then de-spike it.

    Parameters
    ----------
    s, x, y
        Arc length and coordinates of the centerline.
    config
        Curvature settings. ``None`` uses the defaults.

    Returns
    -------
    curvature, theta : ndarray
        Signed curvature (1/m) and per-node tangent angle (rad).
    """
    config = config or CurvatureConfig()
    cur, theta = _CURVATURE_METHODS[config.method](s, x, y)
    if config.despike:
        cur = despike_curvature(cur)
    return cur, theta


def filter_curvature(cur: np.ndarray) -> np.ndarray:
    """Smooth a curvature signal.

    A 5-point second-order Savitzky-Golay pass followed by a 5-point
    first-order pass -- the latter is a plain 5-point moving average.
    """
    if cur.size < 5:
        return np.asarray(cur, dtype=float).copy()
    return savgol_filter(
        savgol_filter(cur, 5, 2, mode="nearest"), 5, 1, mode="nearest")


def phase_lag(cur: np.ndarray, window: int) -> np.ndarray:
    """Impose a downstream phase lag on a curvature signal.

    Replaces the local curvature by a linearly decaying, upstream-only
    weighted average over *window* nodes,

    .. math:: C^{lag}_i = \\sum_{j=0}^{M-1} w_j C_{i-j}, \\qquad
              w_j = \\frac{2}{M} - \\frac{2j}{M(M-1)}

    The weights form a normalised triangular kernel: largest at the local
    node, zero at the upstream end of the window, summing to one. It stands in
    for the exponential memory kernel of linear bend theory, and is what puts
    the pool downstream of the bend apex instead of at it. Near the upstream
    end the window is truncated to the nodes available.

    Parameters
    ----------
    cur
        Curvature signal.
    window
        Window length in nodes. Values below 2 return the input unchanged.

    Returns
    -------
    ndarray
        The phase-lagged curvature.
    """
    cur = np.asarray(cur, dtype=float)
    n = cur.size
    window = int(window)
    if window < 2 or n < 3:
        return cur.copy()

    out = cur.copy()
    # Ramp-up region: the window is shorter than requested, so the kernel
    # changes with position and has to be applied node by node.
    ramp_end = min(window, n)
    for i in range(2, ramp_end):
        m = i
        j = np.arange(m)
        weights = 2.0 / m - j * 2.0 / (m * (m - 1))
        out[i] = float(np.dot(weights, cur[i::-1][:m]))
    # Steady region: one fixed kernel, so the whole thing is a convolution.
    if n > window:
        m = window
        j = np.arange(m)
        weights = 2.0 / m - j * 2.0 / (m * (m - 1))
        full = np.convolve(cur, weights)[:n]
        out[window:] = full[window:]
    return out


# --------------------------------------------------------------------------
# generation and loading
# --------------------------------------------------------------------------

def kinoshita_equation_text(kinoshita: KinoshitaConfig,
                            unicode: bool = True) -> str:
    """Render the Kinoshita equation for the configured parameters."""
    theta0 = kinoshita.theta0
    lam = kinoshita.arc_wavelength
    js, jf = kinoshita.skewness, kinoshita.flatness
    if unicode:
        th, pi = "θ", "π"
    else:
        th, pi = "THETA", "PI"
    text = f"{th} = {theta0:.6g}*sin(2{pi}s/{lam:g})"
    terms = []
    if js:
        terms.append(f"{js:.6g}*cos(6{pi}s/{lam:g})")
    if jf:
        terms.append(f"-{jf:.6g}*sin(6{pi}s/{lam:g})")
    if terms:
        joined = " ".join(terms) if len(terms) == 1 else " ".join(terms)
        text += f" + {theta0 ** 3:.6g}*[{joined}]"
    return text


def build_kinoshita(kinoshita: KinoshitaConfig, channel: ChannelConfig,
                    flip_streamwise: bool = False) -> Centerline:
    """Generate a Kinoshita curve.

    The curve is defined by its direction angle as a function of arc length,

    .. math:: \\theta(s) = \\theta_0 \\sin(2\\pi s/\\lambda)
              + \\theta_0^3 [J_s \\cos(6\\pi s/\\lambda)
              - J_f \\sin(6\\pi s/\\lambda)]

    and the coordinates follow by integrating the tangent.

    Parameters
    ----------
    kinoshita
        Curve parameters.
    channel
        Channel parameters; only ``ds`` is used.
    flip_streamwise
        Reverse the flow direction of the result.

    Returns
    -------
    Centerline
    """
    ds = channel.ds
    lam = kinoshita.arc_wavelength
    theta0 = kinoshita.theta0
    n = int(kinoshita.n_bends * lam / ds) + 1
    s = np.linspace(0.0, kinoshita.n_bends * lam, n)

    theta = (theta0 * np.sin(2 * np.pi * s / lam)
             + theta0 ** 3 * (kinoshita.skewness * np.cos(6 * np.pi * s / lam)
                              - kinoshita.flatness * np.sin(6 * np.pi * s / lam)))
    theta[np.abs(theta) < ZERO] = 0.0

    # x[i] = ds * sum_{j<i} cos(theta[j]); the running sum excludes node i.
    x = ds * np.concatenate(([0.0], np.cumsum(np.cos(theta))[:-1]))
    y = ds * np.concatenate(([0.0], np.cumsum(np.sin(theta))[:-1]))
    x[np.abs(x) < ZERO] = 0.0
    y[np.abs(y) < ZERO] = 0.0

    # One extra node so that the closing segment exists.
    x = np.concatenate((x, [x[-1] + x[1] - x[0]]))
    y = np.concatenate((y, [y[-1] + y[1] - y[0]]))
    s = np.concatenate((s, [s[-1] + ds]))
    theta = np.concatenate((theta, theta[-1:]))

    if flip_streamwise:
        x, y = x[::-1].copy(), y[::-1].copy()
        theta = np.concatenate((theta[::-1][1:], theta[:1]))

    cur = np.zeros(x.size)
    cur[1:n] = np.diff(theta[:n]) / ds
    cur[np.abs(cur) < ZERO] = 0.0
    cur[0], cur[-1] = cur[-2], cur[1]
    return Centerline(x, y, s, cur, theta)


def load_centerline(path: str | Path, flip_streamwise: bool = False
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Read centerline coordinates from a two-column text file.

    Parameters
    ----------
    path
        File with one ``x y`` pair per line, in a projected (metric)
        coordinate system. Blank lines and ``#`` comments are ignored.
    flip_streamwise
        Reverse the node order, i.e. the flow direction.

    Returns
    -------
    x, y : ndarray

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file does not hold at least three coordinate pairs.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"centerline file not found: {path}")
    data = np.loadtxt(path, comments="#")
    data = np.atleast_2d(data)
    if data.shape[1] < 2:
        raise ValueError(
            f"{path}: expected at least two columns (x y), "
            f"found {data.shape[1]}")
    if data.shape[0] < 3:
        raise ValueError(
            f"{path}: expected at least three nodes, found {data.shape[0]}")
    x, y = data[:, 0].astype(float), data[:, 1].astype(float)
    if flip_streamwise:
        x, y = x[::-1].copy(), y[::-1].copy()
    return x, y


def extend_ends(x: np.ndarray, y: np.ndarray, extension: float,
                spacing: float) -> tuple[np.ndarray, np.ndarray]:
    """Add straight, tangent sections at both ends of a centerline.

    Offsetting and the migration model both behave badly at a free end, and a
    flow model fed with the output needs straight inlet and outlet reaches
    anyway.

    The tangent direction is taken as a proper unit vector. v1.x divided the
    end segment by the nominal spacing instead, so its extensions were
    slightly mis-scaled whenever the actual node spacing differed from the
    nominal one -- which, after resampling, it always does a little.

    Parameters
    ----------
    x, y
        Centerline coordinates.
    extension
        Length of each straight section (m).
    spacing
        Node spacing within the straight sections (m).

    Returns
    -------
    x, y : ndarray
        The extended centerline.
    """
    n_ext = max(int(extension / spacing), 1)

    head_len = np.hypot(x[1] - x[0], y[1] - y[0])
    tail_len = np.hypot(x[-1] - x[-2], y[-1] - y[-2])
    if head_len <= 0 or tail_len <= 0:
        log.warning("centerline has a zero-length end segment; "
                    "skipping end extension")
        return np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    hx, hy = (x[1] - x[0]) / head_len, (y[1] - y[0]) / head_len
    tx, ty = (x[-1] - x[-2]) / tail_len, (y[-1] - y[-2]) / tail_len

    back = np.linspace(extension, spacing, n_ext)
    forward = np.linspace(spacing, extension, n_ext)
    head_x, head_y = x[0] - back * hx, y[0] - back * hy
    tail_x, tail_y = x[-1] + forward * tx, y[-1] + forward * ty
    return (np.concatenate((head_x, x, tail_x)),
            np.concatenate((head_y, y, tail_y)))


def build_centerline(config: Config) -> Centerline:
    """Build the initial centerline for *config*, ready for the time loop.

    Dispatches on ``config.mode``, then resamples, smooths, extends the ends
    and computes the curvature.

    Returns
    -------
    Centerline
    """
    channel = config.channel
    if config.mode == "kinoshita":
        line = build_kinoshita(config.kinoshita, channel,
                              config.flip.streamwise)
        x, y = line.x, line.y
        extension, spacing = config.kinoshita.arc_wavelength / 10, channel.ds
    else:
        x, y = load_centerline(config.centerline_file, config.flip.streamwise)
        _, x, y = resample(x, y, channel.interval)
        extension, spacing = channel.width, channel.interval

    x, y = extend_ends(x, y, extension, spacing)
    _, x, y = smooth(x, y, config.curvature.n_passes)
    s, x, y = resample(x, y, channel.interval)
    cur, theta = curvature(s, x, y, config.curvature)
    return Centerline(x, y, s, cur, theta)
