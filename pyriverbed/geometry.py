"""Polyline offsetting: expanding a 1D centerline into a 2D channel.

The bed is computed on a curvilinear ``(s, n)`` grid, but the deliverable is a
point cloud in real ``(x, y)`` coordinates. This module builds it by offsetting
the centerline to either side at every transverse station and attaching the
corresponding bed elevations.
"""

from __future__ import annotations

import numpy as np

from .bed import BedTopography
from .config import ChannelConfig
from .logging_utils import get_logger

log = get_logger(__name__)

__all__ = ["offset_polyline", "build_point_cloud", "bankline_polygon"]

#: Two segments whose normalised dot product exceeds this are treated as
#: collinear, because the intersection they define is numerically useless.
_COLLINEAR = 1.0 - 1e-10


def offset_polyline(x: np.ndarray, y: np.ndarray, distance: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Offset a polyline to both sides by *distance*.

    Displacing each node along its own normal fails at corners: the offset
    self-intersects on the inside of a bend and gaps open on the outside.
    Instead each *segment* is offset as a whole and the offset node is placed
    at the intersection of the two adjacent offset segments, which keeps the
    separation between the two offsets exactly ``2 * distance`` everywhere --
    the constant-width assumption the whole framework rests on.

    Nearly collinear segment pairs would give a singular intersection and fall
    back to the naive normal offset.

    Thanks to Y. Luo for the offsetting method.

    Parameters
    ----------
    x, y
        Centerline coordinates.
    distance
        Offset distance (m).

    Returns
    -------
    xo, yo : ndarray
        Arrays of shape ``(n_nodes, 2)``. Column 0 is the left offset
        (counterclockwise), column 1 the right offset (clockwise).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    dx, dy = np.diff(x), np.diff(y)
    seg = np.hypot(dx, dy)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(seg > 0, distance / seg, 0.0)
    dxl, dyl = dx * scale, dy * scale

    # Naive per-segment offset of the segment's first node, plus a final node
    # carried along the last segment so that every node has a value.
    left_x = np.empty(n)
    left_y = np.empty(n)
    right_x = np.empty(n)
    right_y = np.empty(n)
    left_x[:-1], left_y[:-1] = -dyl + x[:-1], dxl + y[:-1]
    right_x[:-1], right_y[:-1] = dyl + x[:-1], -dxl + y[:-1]
    left_x[-1], left_y[-1] = left_x[-2] + dx[-1], left_y[-2] + dy[-1]
    right_x[-1], right_y[-1] = right_x[-2] + dx[-1], right_y[-2] + dy[-1]

    xo = np.empty((n, 2))
    yo = np.empty((n, 2))
    xo[0], yo[0] = (left_x[0], right_x[0]), (left_y[0], right_y[0])
    xo[-1], yo[-1] = (left_x[-1], right_x[-1]), (left_y[-1], right_y[-1])
    if n < 3:
        return xo, yo

    # Intersect consecutive offset segments. For interior node i the two
    # segments are (i-1) and (i); the 2x2 system has the closed-form solution
    # below, with det = dx[i-1] dy[i] - dx[i] dy[i-1].
    dxa, dya = dx[:-1], dy[:-1]        # segment i-1
    dxb, dyb = dx[1:], dy[1:]          # segment i
    det = dxa * dyb - dxb * dya

    lena, lenb = seg[:-1], seg[1:]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = np.where((lena > 0) & (lenb > 0),
                          (dxa * dxb + dya * dyb) / (lena * lenb), 1.0)
    straight = (cosine > _COLLINEAR) | (np.abs(det) == 0.0)

    for xs, ys, out_x, out_y, column in (
        (left_x, left_y, xo, yo, 0),
        (right_x, right_y, xo, yo, 1),
    ):
        ba = dya * xs[:-2] - dxa * ys[:-2]
        bb = dyb * xs[1:-1] - dxb * ys[1:-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            px = (-ba * dxb + dxa * bb) / det
            py = (dya * bb - dyb * ba) / det
        out_x[1:-1, column] = np.where(straight, xs[1:-1], px)
        out_y[1:-1, column] = np.where(straight, ys[1:-1], py)
    return xo, yo


def build_point_cloud(x: np.ndarray, y: np.ndarray, bed: BedTopography,
                      channel: ChannelConfig,
                      collect_banklines: bool = False):
    """Assemble the riverbed point cloud in ``(x, y, z)`` form.

    Node ordering is centerline first, then the left/right pair at each offset
    station working outwards. The mesh generator relies on this order.

    Parameters
    ----------
    x, y
        Centerline coordinates.
    bed
        Bed topography on the ``(s, n)`` grid.
    channel
        Channel configuration.
    collect_banklines
        Also return the outermost offset pair, i.e. the banklines.

    Returns
    -------
    cloud : ndarray
        Shape ``(n_nodes * (2 * n_offsets + 1), 3)``.
    banklines : tuple of ndarray or None
        ``(left_xy, right_xy)`` when *collect_banklines* is set.
    """
    n_off = channel.n_offsets
    interval = channel.interval
    n = x.size

    blocks = [np.column_stack((x, y, bed.z[:, n_off]))]
    banklines = None
    for i in range(n_off - 1, -1, -1):
        distance = channel.half_width - i * interval
        xo, yo = offset_polyline(x, y, distance)
        blocks.append(np.column_stack((xo[:, 0], yo[:, 0], bed.z[:, -1 - i])))
        blocks.append(np.column_stack((xo[:, 1], yo[:, 1], bed.z[:, i])))
        if i == 0 and collect_banklines:
            banklines = (np.column_stack((xo[:, 0], yo[:, 0])),
                         np.column_stack((xo[:, 1], yo[:, 1])))
    cloud = np.concatenate(blocks, axis=0)
    assert cloud.shape == (n * (2 * n_off + 1), 3)
    return cloud, banklines


def bankline_polygon(banklines: tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    """Close the two banklines into a single polygon.

    Walks up the left bank, back down the right bank, and repeats the first
    node so the ring is closed.
    """
    left, right = banklines
    return np.concatenate((left, right[::-1], left[:1]), axis=0)
