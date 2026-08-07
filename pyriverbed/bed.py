"""Synthetic bed topography for a constant-width meandering channel.

Implements the Beck (1988) analytical bed: the transverse bed slope is taken
proportional to the local curvature, the flow depth then grows linearly
towards the outer bank and decays exponentially towards the inner bank, and
the depth on the centerline follows from conserving the cross-sectional area.

See ``THEORY_GUIDE.md`` section 6 for the derivation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ZERO, ChannelConfig, FlipConfig
from .logging_utils import get_logger

log = get_logger(__name__)

__all__ = ["BedTopography", "transverse_slope", "centerline_depth",
           "compute_bed"]


@dataclass
class BedTopography:
    """Bed elevation on the curvilinear ``(s, n)`` grid.

    Attributes
    ----------
    z
        Bed elevation, shape ``(n_streamwise, 2 * n_offsets + 1)``. Column 0
        is the right bank, column ``n_offsets`` the centerline, the last
        column the left bank. Pools are negative, bars positive.
    transverse_slope
        The transverse bed slope at each section.
    centerline_depth
        The flow depth on the centerline at each section (m).
    """

    z: np.ndarray
    transverse_slope: np.ndarray
    centerline_depth: np.ndarray

    @property
    def n_streamwise(self) -> int:
        """Number of cross sections."""
        return int(self.z.shape[0])

    @property
    def n_transverse(self) -> int:
        """Number of nodes per cross section."""
        return int(self.z.shape[1])

    @property
    def relief(self) -> float:
        """Total bar-to-pool relief (m)."""
        return float(np.max(self.z) - np.min(self.z))


def transverse_slope(cur: np.ndarray, channel: ChannelConfig) -> np.ndarray:
    """Transverse bed slope from the local curvature.

    .. math:: S_T = A H C \\xi_{S_T}

    with the Beck (1988) scour factor ``A``.

    The sign puts the pool against the **outer** bank, which is the physically
    correct arrangement: for a left-turning reach (positive curvature) the
    outer bank is the right-hand one, and that is where the bed is deepest.

    v1.x used the opposite sign here and shipped with its transverse flip
    switched on, which cancelled the error. v2 fixes the sign and defaults the
    flip to off, so the two agree; a v1.x steering file is read with its
    ``FLIPTRANS`` value inverted, which reproduces v1.x output exactly.

    Parameters
    ----------
    cur
        Curvature (1/m).
    channel
        Channel configuration, supplying depth, scour factor and corrector.

    Returns
    -------
    ndarray
        The transverse bed slope at each section.
    """
    return (channel.scour_factor * channel.depth * np.asarray(cur, float)
            * channel.transverse_slope_corrector)


def centerline_depth(st: np.ndarray, channel: ChannelConfig) -> np.ndarray:
    """Flow depth on the centerline, from conservation of cross-sectional area.

    Requiring ``int(zeta) dn = 2 b H`` over the section and solving for the
    centerline depth gives

    .. math:: h_c = \\frac{4 b H |S_T| - S_T^2 b^2}
              {2 b |S_T| + 2H - 2H \\exp(-|S_T| b / H)}

    which tends to ``H`` as ``S_T`` tends to zero, i.e. a flat bed in a
    straight reach. This is what makes the synthetic bed self-consistent:
    deepening a pool automatically raises the opposite bar by the compensating
    amount, so the reach-averaged depth stays put.

    Parameters
    ----------
    st
        Transverse bed slope.
    channel
        Channel configuration.

    Returns
    -------
    ndarray
        Centerline flow depth (m).
    """
    b, h = channel.half_width, channel.depth
    magnitude = np.maximum(np.abs(np.asarray(st, float)), ZERO)
    numerator = 4.0 * b * h * magnitude - magnitude ** 2 * b ** 2
    denominator = (2.0 * b * magnitude + 2.0 * h
                   - 2.0 * h * np.exp(-magnitude * b / h))
    return numerator / denominator


def compute_bed(cur: np.ndarray, s: np.ndarray, channel: ChannelConfig,
                flip: FlipConfig | None = None) -> BedTopography:
    """Build the synthetic bed for a curvature signal.

    Parameters
    ----------
    cur
        Curvature (1/m), one value per cross section. Pass the *phase-lagged*
        curvature to get the bed that a real river has.
    s
        Arc length (m), same length as *cur*, used for the longitudinal fall.
    channel
        Channel configuration.
    flip
        Reflection settings; only ``transverse`` is used.

    Returns
    -------
    BedTopography
    """
    flip = flip or FlipConfig()
    b, h, n_off = channel.half_width, channel.depth, channel.n_offsets

    st = transverse_slope(cur, channel)
    # Keep the slope away from exactly zero so the divisions below stay finite;
    # the resulting profile is flat to within round-off, as it should be.
    st = np.where(np.abs(st) < ZERO, ZERO, st)
    hc = centerline_depth(st, channel)

    n = np.linspace(-b, b, 2 * n_off + 1)          # transverse coordinate
    fall = (np.max(s) - np.asarray(s, float)) * channel.slope

    st_col = st[:, None]
    hc_col = hc[:, None]
    product = st_col * n[None, :]                  # S_T * n

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(n[None, :] != 0.0, hc_col / product, 0.0)
        # Pool side (S_T n < 0): depth grows linearly away from the centerline.
        linear = (1.0 - ratio) * np.maximum(-product, 0.0)
        # Bar side (S_T n > 0): depth decays exponentially.
        exponential = ratio * np.exp(-product / h) * np.maximum(product, 0.0)
    depth = linear + exponential
    depth[:, n_off] = hc                           # the two branches meet here

    z = h - (depth - fall[:, None])
    z[np.abs(z) < ZERO] = 0.0
    if flip.transverse:
        z = z[:, ::-1].copy()
    return BedTopography(z=z, transverse_slope=st, centerline_depth=hc)
