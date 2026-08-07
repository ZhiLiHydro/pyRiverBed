"""Figures and animations.

The three-panel summary figure is the one pyRiverBed has always drawn -- plan
view, curvature signal, curvilinear bed -- because it is what the published
figures show. What is new is that it is built through an explicit Figure/Axes
API instead of the pyplot state machine, so it is safe to call from a GUI
thread or from a notebook, and a diagnostics figure summarising the whole run
has been added.

Series are distinguished by line style as well as by colour, so identity never
rests on colour alone.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import matplotlib
import numpy as np

from .bed import BedTopography
from .config import Config
from .logging_utils import get_logger
from .migration import Cutoff
from .planform import Centerline

log = get_logger(__name__)

__all__ = ["PlotStyle", "RunFrames", "plot_summary", "plot_diagnostics",
           "plot_meander_belt", "plot_centerline_history", "draw_oxbows",
           "add_age_colorbar", "make_gifs"]


@dataclass
class PlotStyle:
    """Colours and sizes of the pyRiverBed figures.

    The plan-view palette deliberately mimics an aerial photograph at low
    flow: submerged bed in dark water-green, exposed bar in pale sand, and
    floodplain in olive. It is what makes a synthetic river recognisable as a
    river.

    Oxbow lakes are shaded by **age** from a single perceptually uniform ramp:
    dark blue-purple for the oldest, yellow for the most recent. Age is the one
    quantity that makes a meander belt readable -- which loops were abandoned
    first, and in what order the belt was built -- and a sequential ramp is the
    right encoding for it. The ends of the ramp are labelled "Old" and "New"
    rather than with step numbers, because the ordering is the message.
    """

    wet: tuple[float, float, float] = (70 / 256, 80 / 256, 60 / 256)
    dry: tuple[float, float, float] = (200 / 256, 187 / 256, 160 / 256)
    background: tuple[float, float, float] = (99 / 256, 120 / 256, 84 / 256)
    oxbow_cmap: str = "plasma"
    oxbow_edge: str = "peru"
    oxbow_edge_width: float = 0.8
    centerline: str = "red"
    belt_centerline: str = "black"
    belt_centerline_width: float = 3.0
    curvature_original: str = "darkorange"
    curvature_filtered: str = "dodgerblue"
    curvature_lagged: str = "orangered"
    bed_cmap: str = "gist_earth"
    figsize: tuple[float, float] = (16.0, 12.0)
    label_size: int = 14
    title_size: int = 16


def _age_norm(cutoffs: Sequence[Cutoff], n_steps: int):
    """Normaliser mapping a cutoff's time step onto the age ramp.

    Spans the whole run rather than only the range of observed cutoffs, so the
    colour of a given oxbow lake does not shift as later ones appear. That is
    what makes the colours comparable between the frames of an animation.
    """
    from matplotlib.colors import Normalize
    top = max(int(n_steps), 1)
    if cutoffs:
        top = max(top, max(cut.step for cut in cutoffs))
    return Normalize(vmin=0, vmax=top)


def draw_oxbows(axes, cutoffs: Sequence[Cutoff], n_steps: int,
                style: PlotStyle | None = None, edge: bool = True,
                zorder: int = 0):
    """Fill the oxbow lakes of *cutoffs*, shaded oldest-to-newest.

    Parameters
    ----------
    axes
        Target axes.
    cutoffs
        The cutoff events to draw.
    n_steps
        Total time steps of the run, used to normalise the age ramp.
    style
        Colours and sizes.
    edge
        Outline each lake, which keeps overlapping lakes legible.
    zorder
        Drawing order.

    Returns
    -------
    matplotlib.cm.ScalarMappable or None
        A mappable suitable for :func:`add_age_colorbar`, or ``None`` if there
        was nothing to draw.
    """
    style = style or PlotStyle()
    if not cutoffs:
        return None
    import matplotlib.cm as cm

    norm = _age_norm(cutoffs, n_steps)
    colours = cm.ScalarMappable(norm=norm, cmap=style.oxbow_cmap)
    # Oldest first, so recent lakes lie on top of the ones they overprint --
    # the same order in which the floodplain was actually built.
    for cut in sorted(cutoffs, key=lambda c: c.step):
        if not cut.oxbow_x.size:
            continue
        colour = colours.to_rgba(cut.step)
        axes.fill(cut.oxbow_x, cut.oxbow_y, color=colour, zorder=zorder)
        if edge:
            axes.plot(cut.oxbow_x, cut.oxbow_y, color=style.oxbow_edge,
                      linewidth=style.oxbow_edge_width, zorder=zorder)
    colours.set_array([])
    return colours


def add_age_colorbar(figure, axes, mappable, label_new: str = "New",
                     label_old: str = "Old", fontsize: int = 12):
    """Attach a thin, tickless "Old to New" colour bar for the oxbow ages.

    Numeric ticks would invite the reader to compare exact ages, which is not
    what the ramp is for; the two words say everything it encodes.
    """
    if mappable is None:
        return None
    bar = figure.colorbar(mappable, ax=axes, fraction=0.035, pad=0.02,
                          aspect=18)
    bar.set_ticks([])
    bar.outline.set_visible(False)
    bar.ax.text(0.5, 1.02, label_new, transform=bar.ax.transAxes,
                ha="center", va="bottom", fontsize=fontsize)
    bar.ax.text(0.5, -0.02, label_old, transform=bar.ax.transAxes,
                ha="center", va="top", fontsize=fontsize)
    return bar


def add_flow_arrow(axes, x: np.ndarray, y: np.ndarray, style: PlotStyle | None
                   = None, label: str = "Flow", fontsize: int = 12):
    """Mark the upstream end with an arrow along the initial flow direction."""
    style = style or PlotStyle()
    if x.size < 2:
        return
    span = max(np.ptp(x), np.ptp(y))
    if span <= 0:
        return
    length = 0.06 * span
    step = max(x.size // 50, 1)
    dx, dy = x[step] - x[0], y[step] - y[0]
    norm = float(np.hypot(dx, dy))
    if norm == 0:
        return
    ux, uy = dx / norm, dy / norm
    tail = (x[0] - length * ux, y[0] - length * uy)
    axes.annotate("", xy=(x[0], y[0]), xytext=tail,
                  arrowprops=dict(arrowstyle="-|>", linewidth=2.2,
                                  color=style.belt_centerline,
                                  shrinkA=0, shrinkB=0))
    # Label above the arrow rather than behind it: the inlet usually sits at the
    # left edge of the data, where a label placed further left lands on top of
    # the y-axis title.
    axes.text(tail[0] + 0.5 * length * ux, tail[1] + 0.5 * length * uy,
              label, ha="center", va="bottom", fontsize=fontsize,
              fontweight="bold", color=style.belt_centerline,
              rotation=np.degrees(np.arctan2(uy, ux)),
              rotation_mode="anchor")


@dataclass
class RunFrames:
    """Bookkeeping of the animation frames written during a run."""

    directories: tuple[Path, Path]
    bed_frames: list[Path] = field(default_factory=list)
    planform_frames: list[Path] = field(default_factory=list)

    def register(self, bed: Path | None, planform: Path | None) -> None:
        """Record one pair of frames."""
        if bed is not None:
            self.bed_frames.append(bed)
        if planform is not None:
            self.planform_frames.append(planform)


def _new_figure(figsize):
    """Create a standalone Figure without touching the pyplot registry."""
    from matplotlib.figure import Figure
    return Figure(figsize=figsize, layout="constrained")


def _save(fig, paths: Sequence[Path], dpi: int) -> None:
    """Render *fig* to every path in *paths*."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    FigureCanvasAgg(fig)
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi)


def plot_summary(
    line: Centerline,
    cloud: np.ndarray,
    bed: BedTopography,
    cur_original: np.ndarray,
    cur_filtered: np.ndarray,
    cur_lagged: np.ndarray,
    config: Config,
    step: int,
    cutoffs: Sequence[Cutoff] = (),
    paths: Sequence[Path] = (),
    dpi: int = 300,
    style: PlotStyle | None = None,
):
    """Draw the three-panel summary figure.

    Panel (a) is the plan view of the synthetic bed, rendered as an exposed-bar
    visualisation at half-bankfull stage, with any oxbow lakes. Panel (b) is
    the curvature signal at its three stages of processing. Panel (c) is the
    bed on the curvilinear grid.

    Parameters
    ----------
    line
        Current centerline.
    cloud
        Riverbed point cloud.
    bed
        Bed topography on the ``(s, n)`` grid.
    cur_original, cur_filtered, cur_lagged
        The curvature signal before filtering, after filtering, and after
        phase lagging.
    config
        Run configuration.
    step
        Current time step.
    cutoffs
        Cutoffs so far, drawn as oxbow lakes.
    paths
        Files to write. Nothing is written if empty.
    dpi
        Output resolution.
    style
        Colours and sizes.

    Returns
    -------
    matplotlib.figure.Figure
    """
    style = style or PlotStyle()
    channel = config.channel
    width, depth, slope = channel.width, channel.depth, channel.slope
    n_off = channel.n_offsets

    fig = _new_figure(style.figsize)
    gs = fig.add_gridspec(8, 8)
    ax_plan = fig.add_subplot(gs[0:6, 0:6])
    ax_cur = fig.add_subplot(gs[6:8, 0:6])
    ax_bed = fig.add_subplot(gs[0:8, 6:8])

    # -- (a) plan view --------------------------------------------------
    title = "(a) Synthetic bed in x-y coordinate\n(exposed sand bar visualization)"
    if config.migration.enabled:
        days = step * config.migration.dt / 86400
        title += f"\nT = {days:.0f} d"
    ax_plan.set_title(title, fontsize=style.title_size)

    draw_oxbows(ax_plan, cutoffs, config.n_steps, style=style, zorder=0)

    water_level = depth / 2
    fall = np.tile((line.s[-1] - line.s) * slope, 2 * n_off + 1)
    elevation = cloud[:, 2] - fall
    wet = elevation < water_level
    ax_plan.scatter(cloud[wet, 0], cloud[wet, 1], color=style.wet, s=2, zorder=1)
    ax_plan.scatter(cloud[~wet, 0], cloud[~wet, 1], color=style.dry, s=2,
                    zorder=1)
    ax_plan.plot(line.x, line.y, color=style.centerline, linewidth=0.5,
                 zorder=2, label="centerline")
    ax_plan.plot(line.x[0], line.y[0], marker="o", color=style.centerline,
                 markersize=8, linestyle="", zorder=3,
                 label="upstream end")
    # Proxy handles so the legend names the two bed states.
    ax_plan.plot([], [], marker="o", color=style.wet, markersize=8,
                 linestyle="", label="river channel")
    ax_plan.plot([], [], marker="o", color=style.dry, markersize=8,
                 linestyle="", label="exposed sand bar")
    # 'datalim' keeps the 1:1 scale while filling the panel, so a long shallow
    # reach does not leave the figure half empty.
    ax_plan.set_aspect("equal", adjustable="datalim")
    ax_plan.set_xlabel("Longitudinal direction (m)", fontsize=style.label_size)
    ax_plan.set_ylabel("Latitudinal direction (m)", fontsize=style.label_size)
    ax_plan.set_facecolor(style.background)
    handles, labels = ax_plan.get_legend_handles_labels()
    order = [labels.index(name) for name in
             ("river channel", "exposed sand bar", "centerline", "upstream end")
             if name in labels]
    ax_plan.legend([handles[i] for i in order], [labels[i] for i in order],
                   edgecolor="black", facecolor="white", framealpha=1,
                   fontsize=style.label_size)

    # -- (b) curvature signal -------------------------------------------
    s_hat = line.s / width
    ax_cur.axhline(0.0, color="0.75", linewidth=0.8, zorder=0)
    ax_cur.plot(s_hat, cur_original * width, ":", color=style.curvature_original,
                linewidth=1, label="original")
    ax_cur.plot(s_hat, cur_filtered * width, "--",
                color=style.curvature_filtered, linewidth=1, label="filtered")
    if config.lag.enabled:
        ax_cur.plot(s_hat, cur_lagged * width, "-",
                    color=style.curvature_lagged, linewidth=1,
                    label="filtered & lagged")
    ax_cur.set_xlim(0, np.max(s_hat) * 1.2)
    if config.migration.enabled:
        ax_cur.set_ylim(-1, 1)
        ax_cur.yaxis.set_major_formatter(
            matplotlib.ticker.FormatStrFormatter("%.1f"))
    ax_cur.legend(edgecolor="black", facecolor="white", framealpha=1)
    ax_cur.set_xlabel(
        "Dimensionless streamwise distance (normalized by channel width)",
        fontsize=style.label_size)
    ax_cur.set_ylabel("Dimensionless curvature\n(normalized by channel width)",
                      fontsize=style.label_size)
    ax_cur.set_title("(b) Curvature signal", fontsize=style.title_size)
    ax_cur.spines[["top", "right"]].set_visible(False)

    # -- (c) bed on the curvilinear grid --------------------------------
    n_row = bed.n_streamwise
    vmax = 1 + slope * np.max(line.s) / depth
    image = ax_bed.imshow(bed.z / depth, cmap=style.bed_cmap,
                          aspect=n_off * 16 / max(n_row, 1),
                          interpolation="bilinear", vmin=-1, vmax=vmax)
    ax_bed.set_xlabel("Transverse direction\n(normalized by channel half-width)",
                      fontsize=style.label_size)
    ax_bed.set_xticks([0, n_off, 2 * n_off], ["-1", "0", "1"])
    ax_bed.set_ylabel("Streamwise direction (normalized by channel length)",
                      fontsize=style.label_size)
    ax_bed.set_yticks(np.linspace(0, n_row, 6),
                      ["0", "0.2", "0.4", "0.6", "0.8", "1"])
    bar = fig.colorbar(image, ax=ax_bed, extend="both", fraction=0.08)
    bar.minorticks_on()
    bar.ax.set_ylabel("Elevation (normalized by water depth)",
                      fontsize=style.label_size)
    ax_bed.set_title("(c) Synthetic bed in s-n coordinate\n"
                     "(bilinearly interpolated)", fontsize=style.title_size)

    if paths:
        _save(fig, [Path(p) for p in paths], dpi)
    return fig


def plot_meander_belt(
    line: Centerline,
    config: Config,
    step: int,
    cutoffs: Sequence[Cutoff] = (),
    paths: Sequence[Path] = (),
    dpi: int = 150,
    style: PlotStyle | None = None,
    history_x: Sequence[np.ndarray] = (),
    history_y: Sequence[np.ndarray] = (),
    show_history: bool = False,
):
    """Draw the meander belt: oxbow lakes shaded by age, current channel on top.

    This is the figure that shows what a migration run actually did. Every
    abandoned loop is filled according to *when* it was abandoned, so the belt
    reads as a stratigraphy: the dark lobes were cut off first, the yellow ones
    most recently, and the black line is where the river is now.

    Parameters
    ----------
    line
        The current centerline.
    config
        Run configuration.
    step
        Current time step.
    cutoffs
        Cutoffs so far, drawn as oxbow lakes shaded by age.
    paths
        Files to write.
    dpi
        Output resolution.
    style
        Colours and sizes.
    history_x, history_y
        Centerline coordinates at each graphic printout, oldest first. Only
        used when *show_history* is set.
    show_history
        Also draw the earlier centerlines as faint grey traces. Off by default,
        because the age-shaded belt already carries the history and the traces
        clutter it.

    Returns
    -------
    matplotlib.figure.Figure
    """
    style = style or PlotStyle()
    fig = _new_figure((14.0, 8.0))
    ax = fig.add_subplot(111)

    mappable = draw_oxbows(ax, cutoffs, config.n_steps, style=style, zorder=0)

    if show_history:
        n_total = max(config.n_steps, 1)
        every = config.migration.plot_every
        for i, (cx, cy) in enumerate(zip(history_x, history_y)):
            shade = min(0.85, 0.35 + 0.5 * (1 - i * every / n_total))
            ax.plot(cx, cy, color=str(shade), linewidth=0.7, zorder=1)

    ax.plot(line.x, line.y, color=style.belt_centerline,
            linewidth=style.belt_centerline_width, solid_capstyle="round",
            zorder=2)
    add_flow_arrow(ax, line.x, line.y, style=style,
                   fontsize=style.label_size)

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("Longitudinal direction (m)", fontsize=style.label_size)
    ax.set_ylabel("Latitudinal direction (m)", fontsize=style.label_size)
    days = step * config.migration.dt / 86400
    ax.set_title(f"{days:,.0f} days", fontsize=style.title_size + 2,
                 fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    add_age_colorbar(fig, ax, mappable, fontsize=style.label_size)

    if paths:
        _save(fig, [Path(p) for p in paths], dpi)
    return fig


#: Kept so that code written against the v2.0 alpha API still works.
plot_centerline_history = plot_meander_belt


def plot_diagnostics(
    sinuosity: np.ndarray,
    migration_rate: np.ndarray,
    cutoffs: Sequence[Cutoff],
    config: Config,
    paths: Sequence[Path] = (),
    dpi: int = 300,
):
    """Draw the time series that tell you whether a migration run behaved.

    Sinuosity and migration rate have entirely different scales, so they get
    one stacked panel each rather than sharing a twinned axis.

    Cutoffs get their own thin panel underneath rather than being ruled across
    the series. Ruling works while cutoffs are rare, but a long run produces
    hundreds of them, and a few hundred vertical lines is an opaque curtain
    that hides the very sawtooth it is supposed to explain. As a rug plot in
    its own lane the timing is still legible at any density, and the two series
    stay clean.

    Parameters
    ----------
    sinuosity, migration_rate
        One value per time step.
    cutoffs
        Cutoff events, drawn as a rug in the bottom panel.
    config
        Run configuration.
    paths
        Files to write.
    dpi
        Output resolution.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from matplotlib.gridspec import GridSpec

    neck = [c.step for c in cutoffs if c.kind == "neck"]
    chute = [c.step for c in cutoffs if c.kind == "chute"]

    fig = _new_figure((11.0, 7.5))
    grid = GridSpec(3, 1, figure=fig, height_ratios=(4.0, 4.0, 1.1),
                    hspace=0.18)
    ax_sin = fig.add_subplot(grid[0])
    ax_rate = fig.add_subplot(grid[1], sharex=ax_sin)
    ax_cut = fig.add_subplot(grid[2], sharex=ax_sin)
    steps = np.arange(sinuosity.size)

    for ax in (ax_sin, ax_rate):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(labelbottom=False)

    ax_sin.plot(steps, sinuosity, color="dodgerblue", linewidth=1.6)
    ax_sin.set_ylabel("Sinuosity")
    ax_sin.set_title("Migration diagnostics", fontsize=15)

    ax_rate.plot(steps, migration_rate, color="seagreen", linewidth=1.6)
    ax_rate.set_ylabel("Mean migration rate\n(m per time step)")

    # The cutoff rug: neck events on the upper row, chute events on the lower,
    # so the two are told apart by position and not by colour. The rows are
    # named by their own tick labels rather than by a legend, which needs no
    # extra space and puts the count where the eye already is.
    rows = ((neck, "0.35", 0.55, "neck"), (chute, "orangered", 0.06, "chute"))
    for events, colour, base, _ in rows:
        if events:
            ax_cut.vlines(events, base, base + 0.39, color=colour,
                          linewidth=0.7, alpha=0.85)
    ax_cut.set_ylim(0, 1)
    ax_cut.set_yticks([base + 0.195 for _, _, base, _ in rows])
    ax_cut.set_yticklabels([f"{name} ({len(events)})"
                            for events, _, _, name in rows], fontsize=9)
    for label, (_, colour, _, _) in zip(ax_cut.get_yticklabels(), rows):
        label.set_color(colour if colour != "0.35" else "0.25")
    ax_cut.tick_params(axis="y", length=0)
    ax_cut.set_title("Cutoffs", fontsize=10, loc="left", pad=2)
    ax_cut.set_xlabel("Time step")
    ax_cut.spines[["top", "right", "left"]].set_visible(False)
    ax_cut.set_xlim(0, max(sinuosity.size - 1, 1))

    if paths:
        _save(fig, [Path(p) for p in paths], dpi)
    return fig


def make_gifs(frames: RunFrames, directory: str | Path, stem: str,
              fps: int = 24) -> list[Path]:
    """Assemble the recorded frames into animated GIFs.

    ``imageio`` is imported here rather than at module level so that the rest
    of pyRiverBed works without it; a missing ``imageio`` costs you the GIFs
    and nothing else.

    Parameters
    ----------
    frames
        The frames recorded during the run.
    directory
        Output directory.
    stem
        Output file stem.
    fps
        Frames per second.

    Returns
    -------
    list of pathlib.Path
        The GIFs written; empty if ``imageio`` is unavailable.
    """
    try:
        import imageio.v3 as iio
    except ImportError:
        try:
            import imageio as iio  # type: ignore[no-redef]
        except ImportError:
            log.warning("imageio is not installed, so no GIF was written. "
                        "Install it with 'pip install imageio', or set "
                        "output.save_gif = no to silence this")
            return []

    directory = Path(directory)
    written: list[Path] = []
    groups = (("0", frames.bed_frames), ("1", frames.planform_frames))
    for suffix, paths in groups:
        if len(paths) < 2:
            continue
        images = [iio.imread(p) for p in paths]
        target = directory / f"{stem}_migration{suffix}.gif"
        try:
            iio.imwrite(target, images, duration=1000 / fps, loop=0)
        except TypeError:  # older imageio signature
            iio.mimsave(target, images, "GIF", fps=fps, loop=0)
        written.append(target)
        log.info("   wrote %s  (%d frames at %d fps)", target.name,
                 len(images), fps)
    return written
