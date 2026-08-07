"""Finite element mesh and boundary condition files.

Because the point cloud built by :mod:`pyriverbed.geometry` is a structured
grid in disguise, the triangulation can be written down directly with no
Delaunay step. Each structured quadrilateral becomes two triangles, and since
the streamwise and transverse node spacings were made equal the triangles come
out close to equilateral.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import ChannelConfig
from .logging_utils import get_logger

log = get_logger(__name__)

__all__ = ["triangulate", "write_xyz", "write_bankline", "write_mesh_files"]


def triangulate(n_row: int, n_offsets: int) -> np.ndarray:
    """Build the triangle connectivity of the structured point cloud.

    Parameters
    ----------
    n_row
        Number of streamwise nodes.
    n_offsets
        Number of polyline offsets per side.

    Returns
    -------
    ndarray
        Shape ``(2 * (n_row - 1) * 2 * n_offsets, 3)`` of 1-based node indices.
    """
    rows = np.arange(n_row - 1)[:, None]           # i
    cols = np.arange(n_offsets)[None, :]           # j
    n = n_row

    inner_a = np.where(cols == 0, rows + 1, rows + 1 + (2 * cols - 1) * n)
    inner_b = np.where(cols == 0, rows + 2, rows + 2 + (2 * cols - 1) * n)
    odd_a = rows + 1 + (2 * cols + 1) * n
    odd_b = rows + 2 + (2 * cols + 1) * n
    even_a = rows + 1 + 2 * cols * n
    even_b = rows + 2 + 2 * cols * n
    next_a = rows + 1 + 2 * (cols + 1) * n
    next_b = rows + 2 + 2 * (cols + 1) * n

    quads = np.stack(
        (
            np.stack((odd_a, inner_a, inner_b), axis=-1),
            np.stack((odd_a, inner_b, odd_b), axis=-1),
            np.stack((even_a, next_a, next_b), axis=-1),
            np.stack((even_a, next_b, even_b), axis=-1),
        ),
        axis=2,
    )
    return quads.reshape(-1, 3)


def _boundary_nodes(n_row: int, n_col: int
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Node indices of the inlet, outlet, left bank and right bank."""
    columns = np.arange(n_col)
    inlet = 1 + columns * n_row
    outlet = (1 + columns) * n_row
    interior = np.arange(1, n_row - 1)
    left = (n_col - 2) * n_row + interior + 1
    right = (n_col - 1) * n_row + interior + 1
    return inlet, outlet, left, right


def write_xyz(cloud: np.ndarray, path: str | Path) -> Path:
    """Write the riverbed point cloud as a 3-column ASCII file.

    Loadable in Blue Kenue, in GIS packages, and by any interpolator.
    """
    path = Path(path)
    np.savetxt(path, cloud, fmt="%.6e")
    log.info("   wrote %s  (%d points)", path.name, cloud.shape[0])
    return path


def write_bankline(polygon: np.ndarray, path: str | Path) -> Path:
    """Write the closed bankline polygon as a Blue Kenue ``.i2s`` file."""
    path = Path(path)
    np.savetxt(path, polygon, fmt="%.6e")
    log.info("   wrote %s  (%d vertices)", path.name, polygon.shape[0])
    return path


def write_mesh_files(cloud: np.ndarray, n_row: int, channel: ChannelConfig,
                     directory: str | Path, stem: str) -> list[Path]:
    """Write the FEM mesh and boundary condition files.

    Four files are produced:

    ``<stem>_mesh.t3s``
        Blue Kenue T3 mesh, the input for building a TELEMAC Selafin geometry.
    ``<stem>_mesh.dat``
        Tecplot ``FETRIANGLE`` zone, for visualisation.
    ``<stem>_BC.cli``
        TELEMAC boundary conditions, usable directly.
    ``<stem>_BC.bc2``
        Blue Kenue boundary conditions, for inspecting and editing BC codes.

    The inlet is tagged as prescribed-discharge with a free surface and the
    outlet as prescribed-elevation; the banks stay closed walls. Those are the
    conventional choices for a steady subcritical run, and a starting point
    rather than a prescription.

    Parameters
    ----------
    cloud
        The point cloud, ordered as :func:`~pyriverbed.geometry.build_point_cloud`
        produces it.
    n_row
        Number of streamwise nodes.
    channel
        Channel configuration, for the number of offsets.
    directory
        Output directory.
    stem
        Output file stem.

    Returns
    -------
    list of pathlib.Path
        The files written.
    """
    directory = Path(directory)
    n_col = 2 * channel.n_offsets + 1
    n_node = cloud.shape[0]
    triangles = triangulate(n_row, channel.n_offsets)
    n_ele = triangles.shape[0]
    written: list[Path] = []

    t3s = directory / f"{stem}_mesh.t3s"
    with open(t3s, "w", encoding="utf-8") as f:
        f.write(f":NodeCount {n_node}\n:ElementCount {n_ele}\n#\n:EndHeader\n")
        np.savetxt(f, cloud, fmt="%.6e")
        np.savetxt(f, triangles, fmt="%d")
        f.write("\n\n")
    written.append(t3s)

    dat = directory / f"{stem}_mesh.dat"
    with open(dat, "w", encoding="utf-8") as f:
        f.write(f'TITLE = "{stem}_mesh"\n'
                f'VARIABLES = "X", "Y", "{stem}_mesh"\n'
                f"ZONE NODES={n_node}, ELEMENTS={n_ele}, "
                f"DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n")
        np.savetxt(f, cloud, fmt="%.6e")
        np.savetxt(f, triangles, fmt="%d")
        f.write("\n\n")
    written.append(dat)

    inlet, outlet, left, right = _boundary_nodes(n_row, n_col)
    boundary = np.concatenate((inlet, outlet, left, right))
    n_bnd = boundary.size
    cli = np.zeros((n_bnd, 13), dtype=int)
    cli[:, :2] = 2
    cli[:, 7] = 2
    cli[:, 11] = boundary
    cli[:, 12] = np.arange(n_bnd) + 1
    cli[:n_col, 0], cli[:n_col, 1], cli[:n_col, 2] = 4, 5, 5
    cli[:n_col, 7] = 4
    cli[n_col:2 * n_col, 0] = 5
    cli[n_col:2 * n_col, 1], cli[n_col:2 * n_col, 2] = 4, 4
    cli[n_col:2 * n_col, 7] = 4

    cli_path = directory / f"{stem}_BC.cli"
    with open(cli_path, "w", encoding="utf-8") as f:
        for i, row in enumerate(cli):
            if i < n_col:
                tag = " #Inlet"
            elif i < 2 * n_col:
                tag = " #Outlet"
            else:
                tag = " #"
            f.write(" ".join(str(v) for v in row) + tag + "\n")
        f.write("\n")
    written.append(cli_path)

    bc2 = directory / f"{stem}_BC.bc2"
    with open(bc2, "w", encoding="utf-8") as f:
        f.write(
            ":FileType bc2  ASCII  EnSim 1.0"
            f"\n:NodeCount {n_node}"
            f"\n:ElementCount {n_ele}"
            "\n:ElementType  T3"
            "\n:BoundarySegmentCount 2"
            "\n# id  code  sectionCount startNode1 endNode1 startNode2"
            " endNode2 tracerCode name"
            f'\n:BoundarySegment 1  455  1 1 {n_col} 1 1  4  "Inlet"'
            f'\n:BoundarySegment 2  544  1 {n_col + 1} {2 * n_col} 1 1  4'
            ' "Outlet"'
            "\n:ShorelineCount 1"
            f"\n:ShorelineNodeCount {n_bnd}"
            "\n:EndHeader"
            f"\n:BeginNodes {n_node}\n")
        flat = cloud.copy()
        flat[:, 2] = 0.0
        np.savetxt(f, flat, fmt="%.6e")
        f.write(f":EndNodes\n:BeginElements {n_ele}\n")
        np.savetxt(f, triangles, fmt="%d")
        f.write(f":EndElements\n:BeginTable {n_bnd} 15\n")
        f.write(cli_path.read_text(encoding="utf-8")[:-1])
        f.write(":EndTable\n\n")
    written.append(bc2)

    log.info("   wrote %s  (%d nodes, %d triangles)",
             ", ".join(p.name for p in written), n_node, n_ele)
    return written
