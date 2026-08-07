"""Jupyter notebook frontend for pyRiverBed.

A third way to drive the same model, alongside the CLI and the GUI. Everything
here is a convenience wrapper: nothing in this module is required to use
pyRiverBed from a notebook, it just removes the boilerplate.

Import it as::

    from pyriverbed.notebook import quick_run, show, configure

``ipywidgets`` is used **if it happens to be installed** -- :func:`form` then
gives you an interactive parameter panel -- but it is not a dependency and
everything else works without it.

Examples
--------
>>> from pyriverbed.notebook import quick_run, show      # doctest: +SKIP
>>> result = quick_run(mode='kinoshita', n_bends=4, width=0.6, depth=0.15)
>>> show(result)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

from ._version import __version__
from .config import Config, ConfigError, default_config, dump_ini
from .logging_utils import get_logger, setup_logging
from .model import RiverBedModel, RunResult

log = get_logger(__name__)

__all__ = ["configure", "quick_run", "show", "show_planform", "show_bed",
           "show_curvature", "show_diagnostics", "show_art", "art_styles",
           "save_art", "cross_section", "form", "print_config", "SHORTCUTS"]


#: Keyword shortcuts accepted by :func:`configure` and :func:`quick_run`,
#: mapping a flat name onto its ``section.key`` path. They exist so that a
#: notebook cell reads ``width=0.6`` instead of ``config.channel.width = 0.6``.
SHORTCUTS: dict[str, str] = {
    "mode": "mode",
    "centerline_file": "centerline_file",
    "name": "name",
    # channel
    "width": "channel.width",
    "depth": "channel.depth",
    "slope": "channel.slope",
    "st_corrector": "channel.transverse_slope_corrector",
    "ds": "channel.ds",
    "n_offsets": "channel.n_offsets",
    # kinoshita
    "n_bends": "kinoshita.n_bends",
    "arc_wavelength": "kinoshita.arc_wavelength",
    "theta0": "kinoshita.max_angular_amplitude",
    "skewness": "kinoshita.skewness",
    "flatness": "kinoshita.flatness",
    # curvature and lag
    "curvature_method": "curvature.method",
    "smoothing": "curvature.smoothing_level",
    "migration_smoothing": "curvature.migration_smoothing_level",
    "lag": "lag.enabled",
    "lag_strength": "lag.strength",
    # flip
    "flip_streamwise": "flip.streamwise",
    "flip_transverse": "flip.transverse",
    # migration
    "migration": "migration.enabled",
    "n_steps": "migration.n_steps",
    "dt": "migration.dt",
    "e0": "migration.e0",
    "ub0": "migration.ub0",
    "c0": "migration.c0",
    "cf0": "migration.cf0",
    "fr0": "migration.fr0",
    "end_taper": "migration.end_taper_widths",
    "seed": "migration.seed",
    "log_every": "migration.log_every",
    "plot_every": "migration.plot_every",
    # cutoffs
    "neck_cutoff": "neck_cutoff.enabled",
    "neck_end_margin": "neck_cutoff.end_margin_widths",
    "chute_cutoff": "chute_cutoff.enabled",
    "chute_frequency": "chute_cutoff.frequency",
    "chute_start": "chute_cutoff.start_step",
    "chute_entrance": "chute_cutoff.entrance",
    "chute_span": "chute_cutoff.span",
    "chute_max_angle": "chute_cutoff.max_valley_angle",
    "chute_min_length": "chute_cutoff.min_length_widths",
    "chute_min_sinuosity": "chute_cutoff.min_sinuosity",
    # output
    "output_dir": "output.directory",
    "save_xyz": "output.save_xyz",
    "save_bankline": "output.save_bankline",
    "save_mesh": "output.save_mesh",
    "save_figures": "output.save_figures",
    "save_gif": "output.save_gif",
    "gif_fps": "output.gif_fps",
    "log_file": "output.log_file",
}


def _apply(config: Config, key: str, value: Any) -> None:
    """Set one ``section.key`` path on *config*."""
    path = SHORTCUTS.get(key, key)
    parts = path.split(".")
    target = config
    for part in parts[:-1]:
        target = getattr(target, part)
    if not hasattr(target, parts[-1]):
        raise ConfigError(
            f"unknown parameter {key!r}; the available shortcuts are "
            f"{', '.join(sorted(SHORTCUTS))}")
    setattr(target, parts[-1], value)


def configure(base: Config | str | Path | None = None, **overrides) -> Config:
    """Build a configuration from flat keyword arguments.

    Parameters
    ----------
    base
        Configuration or input file to start from. ``None`` uses the defaults.
    **overrides
        Any name in :data:`SHORTCUTS`, or a full dotted path such as
        ``chute_cutoff__frequency`` written as ``'chute_cutoff.frequency'``.

    Returns
    -------
    Config
        The validated configuration.

    Examples
    --------
    >>> config = configure(width=0.8, depth=0.2, n_bends=5)   # doctest: +SKIP
    """
    if base is None:
        config = default_config()
    elif isinstance(base, Config):
        config = base
    else:
        from .config import load_config
        config = load_config(base)
    for key, value in overrides.items():
        _apply(config, key, value)
    return config.validate()


def quick_run(base: Config | str | Path | None = None,
              verbose: bool = True, **overrides) -> RunResult:
    """Configure and run pyRiverBed in one call.

    Output goes to the current directory unless ``output_dir`` says otherwise,
    and file writing can be switched off entirely to keep a notebook clean::

        result = quick_run(n_bends=4, save_mesh=False, save_xyz=False)

    Parameters
    ----------
    base
        Configuration or input file to start from.
    verbose
        Show the run log. ``False`` shows warnings and errors only.
    **overrides
        As for :func:`configure`.

    Returns
    -------
    RunResult
    """
    import sys
    setup_logging(level=logging.INFO if verbose else logging.WARNING,
                  colour=False, stream=sys.stdout)
    config = configure(base, **overrides)
    return RiverBedModel(config).run()


def print_config(config: Config) -> None:
    """Print a configuration in the input file format."""
    print(dump_ini(config, header=False))


# --------------------------------------------------------------------------
# display helpers
# --------------------------------------------------------------------------

def _pyplot():
    """Import pyplot lazily, so importing this module stays cheap."""
    import matplotlib.pyplot as plt
    return plt


def show(result: RunResult, figsize: tuple[float, float] = (14, 10)):
    """Show the full three-panel summary figure for a finished run.

    Parameters
    ----------
    result
        The run to display.
    figsize
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from . import plotting
    plt = _pyplot()
    style = plotting.PlotStyle(figsize=figsize)
    figure = plotting.plot_summary(
        result.centerline, result.cloud, result.bed,
        result.curvature_original, result.curvature_filtered,
        result.curvature_lagged, result.config, result.steps_completed,
        cutoffs=result.cutoffs, style=style)
    plt.close(figure)
    return figure


def show_planform(result: RunResult, figsize: tuple[float, float] = (12, 7),
                  show_oxbows: bool = True, show_initial: bool = True):
    """Plot the meander belt: oxbow lakes shaded by age, current channel on top.

    Oxbow lakes are filled according to *when* they were abandoned -- dark for
    the oldest, yellow for the most recent -- so the belt reads as a
    stratigraphy rather than an undifferentiated blob.

    Parameters
    ----------
    result
        The run to display.
    figsize
        Figure size in inches.
    show_oxbows
        Fill in the oxbow lakes left by cutoffs.
    show_initial
        Also draw the initial centerline as a dashed grey line.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from . import plotting
    plt = _pyplot()
    style = plotting.PlotStyle()
    figure, axes = plt.subplots(figsize=figsize, layout="constrained")

    mappable = None
    if show_oxbows:
        mappable = plotting.draw_oxbows(axes, result.cutoffs,
                                       result.config.n_steps, style=style,
                                       zorder=0)
    if show_initial and result.centerline_history:
        first_x, first_y = result.centerline_history[0]
        axes.plot(first_x, first_y, color="0.55", linewidth=1, linestyle="--",
                  label="initial centerline", zorder=1)
    line = result.centerline
    axes.plot(line.x, line.y, color=style.belt_centerline, linewidth=2.2,
              solid_capstyle="round", label="current centerline", zorder=2)
    plotting.add_flow_arrow(axes, line.x, line.y, style=style, fontsize=11)

    axes.set_aspect("equal", adjustable="datalim")
    axes.set_xlabel("Longitudinal direction (m)")
    axes.set_ylabel("Latitudinal direction (m)")
    axes.set_title(f"Planform after {result.steps_completed} steps  |  "
                   f"sinuosity {line.sinuosity:.3f}  |  "
                   f"{len(result.cutoffs)} cutoff(s)")
    axes.legend(frameon=True, edgecolor="black", facecolor="white",
                framealpha=1, loc="upper right", fontsize=9)
    axes.spines[["top", "right"]].set_visible(False)
    plotting.add_age_colorbar(figure, axes, mappable, fontsize=11)
    return figure


def show_bed(result: RunResult, figsize: tuple[float, float] = (6, 9)):
    """Plot the bed on the curvilinear ``(s, n)`` grid.

    Returns
    -------
    matplotlib.figure.Figure
    """
    plt = _pyplot()
    config = result.config
    depth = config.channel.depth
    n_off = config.channel.n_offsets
    figure, axes = plt.subplots(figsize=figsize, layout="constrained")
    image = axes.imshow(result.bed.z / depth, cmap="gist_earth",
                        aspect=n_off * 16 / max(result.bed.n_streamwise, 1),
                        interpolation="bilinear", vmin=-1,
                        vmax=1 + config.channel.slope
                        * np.max(result.centerline.s) / depth)
    axes.set_xticks([0, n_off, 2 * n_off], ["-1", "0", "1"])
    axes.set_xlabel("Transverse direction\n(normalized by channel half-width)")
    axes.set_ylabel("Streamwise node")
    bar = figure.colorbar(image, ax=axes, extend="both", fraction=0.08)
    bar.ax.set_ylabel("Elevation (normalized by water depth)")
    axes.set_title("Synthetic bed in s-n coordinate")
    return figure


def show_curvature(result: RunResult, figsize: tuple[float, float] = (12, 4)):
    """Plot the curvature signal at its three stages of processing.

    Returns
    -------
    matplotlib.figure.Figure
    """
    plt = _pyplot()
    width = result.config.channel.width
    s_hat = result.centerline.s / width
    figure, axes = plt.subplots(figsize=figsize, layout="constrained")
    axes.axhline(0, color="0.8", linewidth=0.8)
    axes.plot(s_hat, result.curvature_original * width, ":",
              color="darkorange", linewidth=1, label="original")
    axes.plot(s_hat, result.curvature_filtered * width, "--",
              color="dodgerblue", linewidth=1, label="filtered")
    if result.config.lag.enabled:
        axes.plot(s_hat, result.curvature_lagged * width, "-",
                  color="orangered", linewidth=1, label="filtered & lagged")
    axes.set_xlabel("Dimensionless streamwise distance (channel widths)")
    axes.set_ylabel("Dimensionless curvature")
    axes.set_title("Curvature signal")
    axes.legend(frameon=True, edgecolor="black", facecolor="white",
                framealpha=1)
    axes.spines[["top", "right"]].set_visible(False)
    return figure


def show_diagnostics(result: RunResult):
    """Plot the sinuosity and migration rate time series of a migration run.

    Returns
    -------
    matplotlib.figure.Figure or None
        ``None`` if the run had no time stepping.
    """
    if result.sinuosity.size == 0:
        log.warning("this run had no migration, so there is nothing to plot")
        return None
    from . import plotting
    plt = _pyplot()
    figure = plotting.plot_diagnostics(result.sinuosity, result.migration_rate,
                                       result.cutoffs, result.config)
    plt.close(figure)
    return figure


def show_art(result: RunResult, style_name: str = "fisk", **kwargs):
    """Render one art print of *result* and return it for inline display.

    Art prints have no axes and no legend: the river is a filled ribbon on a
    flat ground, and the composition is all there is. ``'fisk'`` is a homage to
    Harold Fisk's 1944 maps of the Mississippi meander belt and needs a
    migration run, since it draws every historical course of the river.

    Parameters
    ----------
    result
        A finished run.
    style_name
        One of :data:`pyriverbed.art.STYLES`. Call :func:`art_styles` to list
        them.
    **kwargs
        Passed to :class:`~pyriverbed.art.ArtStyle`, so ``paper='a3'``,
        ``dpi=150``, ``title='My River'`` and the rest all work here.

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> show_art(result, "fisk", dpi=120)              # doctest: +SKIP
    >>> show_art(result, "nocturne", title="")         # doctest: +SKIP
    """
    from . import art as art_module
    kwargs.setdefault("dpi", 120)      # screen, not print
    plt = _pyplot()
    figure = art_module.render(result, style_name,
                               art_module.ArtStyle(**kwargs))
    plt.close(figure)
    return figure


def art_styles() -> dict[str, str]:
    """The art styles and a one-line description of each."""
    from . import art as art_module
    return {name: (function.__doc__ or "").strip().splitlines()[0]
            for name, function in sorted(art_module.STYLES.items())}


def save_art(result: RunResult, directory: str | Path = ".",
             prefix: str = "art", styles: Any = None, **kwargs):
    """Render every art style of *result* into *directory*.

    Parameters
    ----------
    result
        A finished run.
    directory
        Where to write. Created if missing.
    prefix
        File name stem; the style name is appended.
    styles
        Style names to render. ``None`` renders all of them.
    **kwargs
        Passed to :class:`~pyriverbed.art.ArtStyle`.

    Returns
    -------
    list of Path
        The files written.
    """
    from . import art as art_module
    return art_module.save_gallery(result, directory, prefix=prefix,
                                   styles=styles,
                                   style=art_module.ArtStyle(**kwargs))


def cross_section(result: RunResult, fraction: float = 0.5,
                  figsize: tuple[float, float] = (8, 4)):
    """Plot one cross section of the bed, at *fraction* along the reach.

    Shows the asymmetry that is the whole point of the Beck profile: linear on
    the pool side, exponential on the bar side.

    Parameters
    ----------
    result
        The run to display.
    fraction
        Position along the reach, from 0 (upstream) to 1 (downstream).
    figsize
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    plt = _pyplot()
    config = result.config
    channel = config.channel
    index = int(np.clip(fraction, 0, 1) * (result.bed.n_streamwise - 1))
    n = np.linspace(-channel.half_width, channel.half_width,
                    2 * channel.n_offsets + 1)
    z = result.bed.z[index]

    # z = 0 is the bed of an equivalent straight channel, so the water surface
    # sits one reach-averaged depth above it, plus the longitudinal fall.
    fall = (result.centerline.s[-1] - result.centerline.s[index]) * channel.slope
    surface = channel.depth + fall
    floor = min(z.min(), 0.0) - 0.15 * (z.max() - z.min() + channel.depth)

    figure, axes = plt.subplots(figsize=figsize, layout="constrained")
    axes.fill_between(n, floor, z, color="#c8bba0", zorder=1, label="sediment")
    axes.fill_between(n, z, surface, where=z < surface, color="#9ecae1",
                      alpha=0.85, zorder=2, label="water")
    axes.plot(n, z, color="#5b4636", linewidth=2, zorder=3)
    axes.axhline(surface, color="#2471a3", linewidth=1.2, zorder=4,
                 label="water surface")
    axes.axhline(fall, color="0.45", linewidth=1.0, linestyle="--", zorder=4,
                 label="straight-channel bed level")
    axes.set_xlim(n[0], n[-1])
    axes.set_ylim(floor, surface + 0.1 * (surface - floor))
    axes.set_xlabel("Transverse coordinate n (m)   "
                    "(negative = right bank, positive = left bank)")
    axes.set_ylabel("Elevation (m)")
    axes.set_title(f"Cross section at {fraction:.0%} of the reach  |  "
                   f"dimensionless curvature "
                   f"{result.curvature_lagged[index] * channel.width:+.3f}")
    axes.legend(frameon=True, edgecolor="black", facecolor="white",
                framealpha=1, fontsize=9, loc="lower left")
    axes.spines[["top", "right"]].set_visible(False)
    return figure


# --------------------------------------------------------------------------
# optional interactive form
# --------------------------------------------------------------------------

def form(config: Config | RunResult | None = None):
    """Return an ``ipywidgets`` parameter panel with a Run button.

    Requires ``ipywidgets``, which is *not* a pyRiverBed dependency. Without
    it, a message explains how to install it and the rest of the notebook API
    remains available.

    Parameters
    ----------
    config
        Configuration to pre-load. A :class:`~pyriverbed.model.RunResult` is
        also accepted, so a form can be seeded from a previous run.

    Returns
    -------
    ipywidgets.Widget or None
        The panel, or ``None`` if ``ipywidgets`` is missing. The most recent
        result is left on the returned object as ``panel.result``.
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError:
        print("ipywidgets is not installed, so the interactive form is "
              "unavailable.\nInstall it with 'pip install ipywidgets', or use "
              "quick_run(...) instead:\n\n"
              "    from pyriverbed.notebook import quick_run, show\n"
              "    result = quick_run(n_bends=4, width=0.6, depth=0.15)\n"
              "    show(result)")
        return None

    if isinstance(config, RunResult):
        config = config.config
    config = config or default_config()
    fields = {
        "mode": widgets.Dropdown(options=list(Config.MODES), value=config.mode,
                                 description="mode"),
        "centerline_file": widgets.Text(value=config.centerline_file,
                                        description="file"),
        "width": widgets.FloatText(value=config.channel.width,
                                   description="width (m)"),
        "depth": widgets.FloatText(value=config.channel.depth,
                                   description="depth (m)"),
        "n_bends": widgets.IntSlider(value=config.kinoshita.n_bends, min=1,
                                     max=20, description="bends"),
        "smoothing": widgets.IntSlider(value=config.curvature.smoothing_level,
                                       min=0, max=100, description="smoothing"),
        "lag_strength": widgets.FloatSlider(value=config.lag.strength, min=0.5,
                                            max=12, step=0.5,
                                            description="lag (widths)"),
        "migration": widgets.Checkbox(value=config.migration.enabled,
                                      description="migration"),
        "n_steps": widgets.IntText(value=config.migration.n_steps,
                                   description="steps"),
        "chute_cutoff": widgets.Checkbox(value=config.chute_cutoff.enabled,
                                         description="chute cutoffs"),
        "chute_frequency": widgets.FloatSlider(
            value=config.chute_cutoff.frequency, min=0.0, max=1.0, step=0.01,
            description="chute freq"),
        "output_dir": widgets.Text(value=config.output.directory,
                                   description="output"),
    }
    run_button = widgets.Button(description="Run pyRiverBed",
                                button_style="primary", icon="play")
    output = widgets.Output()
    panel = widgets.VBox([
        widgets.HTML(f"<b>pyRiverBed {__version__}</b>"),
        widgets.GridBox(list(fields.values()),
                        layout=widgets.Layout(
                            grid_template_columns="repeat(2, 360px)")),
        run_button,
        output,
    ])
    panel.result = None

    def on_click(_):
        output.clear_output()
        with output:
            values = {key: widget.value for key, widget in fields.items()}
            try:
                panel.result = quick_run(**values)
            except ConfigError as exc:
                print(f"invalid input: {exc}")
                return
            display(show(panel.result))

    run_button.on_click(on_click)
    display(panel)
    return panel
