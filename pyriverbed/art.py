"""Art prints from a pyRiverBed run.

:mod:`pyriverbed.plotting` draws *figures*: labelled axes, colour bars,
legends, everything a reader needs to check a number. This module draws
*prints*. No axes, no ticks, no legend -- the channel is a filled ribbon on a
flat ground, and the composition is the whole content.

The flagship style is a homage to Harold Fisk's 1944 maps of the Mississippi
River meander belt, made for the US Army Corps of Engineers and now the most
reproduced piece of fluvial geomorphology there is. Fisk's insight was to draw
every historical course of the river in its own flat colour and let them
overlap: the 1880 course in green, 1820 in salmon-pink, 1765 in light blue, and
behind them the prehistoric courses, a palimpsest of meanders on cream paper.

That is exactly the data a migration run produces. ``centerline_history`` is a
sequence of channel courses and every :class:`~pyriverbed.migration.Cutoff`
carries an oxbow lake stamped with the step it was abandoned, so a run can be
drawn the way Fisk drew a survey.

Styles
------
``fisk``
    Overlapping flat-colour historical courses on aged paper, with a title
    cartouche and a scale bar. The Fisk homage.
``strata``
    The same courses in one perceptual ramp, oldest to newest, on a dark
    ground. Reads as a stratigraphic section seen from above.
``blueprint``
    Cyanotype: pale channel on Prussian blue, hairline grid, drafting
    annotations.
``nocturne``
    Warm metallic channel with a soft glow on near-black. The gold-foil look.
``minimal``
    One accent ribbon on a large flat ground, generous margins, tiny caption.
``bathymetry``
    The bed itself, as filled depth bands with crisp banklines.
``contour``
    The bed as thin contour lines only. Topographic line art.

Usage
-----
>>> from pyriverbed import art                          # doctest: +SKIP
>>> art.render(result, "fisk")                          # doctest: +SKIP
>>> art.save_gallery(result, "img", prefix="v2_art")    # doctest: +SKIP

Nothing here needs anything beyond NumPy and Matplotlib.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

from .geometry import bankline_polygon, offset_polyline
from .logging_utils import get_logger

log = get_logger(__name__)

__all__ = ["ArtStyle", "STYLES", "PAPERS", "render", "save_gallery",
           "channel_ribbon", "auto_window", "fisk", "strata", "blueprint",
           "nocturne", "minimal", "bathymetry", "contour"]


# --------------------------------------------------------------------------
# palettes and paper sizes
# --------------------------------------------------------------------------

#: Flat, saturated course colours in the spirit of Fisk's sheets. Ordered
#: oldest to newest, so a run with more courses than colours cycles through
#: them the way Fisk reused inks across a sheet.
FISK_COURSES: tuple[str, ...] = (
    "#7C6A9C",   # plum
    "#4E7A3E",   # forest green
    "#C9662E",   # burnt orange
    "#5B8CA8",   # dusty blue
    "#D8A94B",   # ochre
    "#8C4A52",   # brick
    "#3F7A72",   # teal
    "#B5763F",   # tan
    "#6B7F3A",   # olive
    "#A85B7A",   # mulberry
)

#: Aged-paper ground for the Fisk style.
FISK_PAPER = "#EFE3C8"
FISK_INK = "#3A3128"

#: Hatch patterns, cycled across the older courses. Fisk overprinted his
#: colours with stipple and hatching so that overlapping courses stayed
#: legible where their colours were similar; the same trick keeps a run with
#: many courses readable, and it survives being printed in greyscale.
FISK_HATCHES: tuple[str, ...] = ("", "....", "", "///", "", "\\\\\\", "",
                                 "...", "", "//")

#: Print sizes in inches, portrait. ``render`` transposes them when the reach
#: is wider than it is tall, which meander belts almost always are.
#:
#: ``fit`` is not a paper size but an instruction: shape the canvas to the reach
#: itself. A meander belt is typically four or five times wider than it is tall,
#: so on any standard sheet it becomes a thin band across the middle with half
#: the paper left empty. ``fit`` is the default for that reason -- panoramic is
#: the honest format for this subject -- and the named sizes are there for when
#: the print has to go in a standard frame.
PAPERS: dict[str, tuple[float, float]] = {
    "a4": (8.27, 11.69),
    "a3": (11.69, 16.54),
    "a2": (16.54, 23.39),
    "letter": (8.5, 11.0),
    "tabloid": (11.0, 17.0),
    "square": (12.0, 12.0),
}

#: Long edge, in inches, of a ``paper = "fit"`` canvas.
FIT_LONG_EDGE = 18.0

#: Aspect ratios beyond these are clamped, so a nearly straight reach does not
#: produce a canvas too extreme to print or to look at.
FIT_ASPECT_LIMITS = (0.35, 3.2)


@dataclass
class ArtStyle:
    """Knobs shared by every print.

    Attributes
    ----------
    paper
        ``'fit'`` to shape the canvas to the reach, a key into :data:`PAPERS`
        for a standard sheet, or an explicit ``(width, height)`` in inches.
    dpi
        Output resolution. 300 is print quality; 150 is plenty on screen.
    margin
        Fraction of the shorter side left blank around the artwork.
    n_courses
        How many historical courses to draw. The history is subsampled evenly
        if it holds more than this. Beyond a dozen or so the overlaps stop
        being readable.
    show_oxbows
        Draw the abandoned loops as lakes.
    show_current
        Draw the final course on top of the historical ones.
    title, subtitle
        Cartouche text. ``None`` derives them from the run; ``""`` omits the
        cartouche entirely.
    scale_bar
        Draw a scale bar.
    texture
        Strength of the paper mottling, 0 to 1. 0 switches it off.
    bed_grid
        Nominal grid resolution used by the ``bathymetry`` and ``contour``
        styles, which have to resample the point cloud before they can contour
        it. Cost grows as the square, so this is the knob to turn down when
        iterating on a composition and up for the final print.
    window
        ``(start, end)`` fractions of the reach to show, or ``None`` for all of
        it.

        The belt styles want the whole reach. The two bed styles do not: after
        a long migration run a reach is hundreds of channel widths long, and at
        that zoom the channel is a hairline and the bar and pool the print is
        supposed to show are invisible. They therefore choose a window of about
        45 channel widths for themselves unless this says otherwise.
    seed
        Seed for the texture, so a print is reproducible.
    """

    paper: str | tuple[float, float] = "fit"
    dpi: int = 300
    margin: float = 0.06
    n_courses: int = 8
    show_oxbows: bool = True
    show_current: bool = True
    title: str | None = None
    subtitle: str | None = None
    scale_bar: bool = True
    texture: float = 0.35
    bed_grid: int = 700
    window: tuple[float, float] | None = None
    seed: int | None = 0

    def figsize(self, aspect: float) -> tuple[float, float]:
        """Canvas size in inches for a reach of width-to-height ratio *aspect*.

        For ``paper = 'fit'`` the canvas takes the reach's own aspect ratio, so
        the artwork fills the frame instead of sitting in a band across a sheet
        it does not match. Otherwise a named or explicit size is used, turned
        landscape when the reach is wider than it is tall.
        """
        if self.paper == "fit":
            low, high = FIT_ASPECT_LIMITS
            aspect = float(np.clip(aspect, low, high))
            return ((FIT_LONG_EDGE, FIT_LONG_EDGE / aspect) if aspect >= 1.0
                    else (FIT_LONG_EDGE * aspect, FIT_LONG_EDGE))
        size = (PAPERS[self.paper] if isinstance(self.paper, str)
                else tuple(self.paper))
        short, long = min(size), max(size)
        return (long, short) if aspect >= 1.0 else (short, long)


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------

def channel_ribbon(x: np.ndarray, y: np.ndarray, half_width: float
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Closed polygon of a channel of half-width *half_width* about a course.

    A river drawn as a line is a graph; drawn as a ribbon of its own width it
    reads as a river. Uses the same segment-intersection offsetting as the
    model itself, so the ribbon is exactly the channel the model computed.

    Parameters
    ----------
    x, y
        Centerline coordinates.
    half_width
        Half the channel width (m).

    Returns
    -------
    px, py : ndarray
        Closed polygon, first node repeated at the end.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3:
        return x.copy(), y.copy()
    xo, yo = offset_polyline(x, y, half_width)
    left = np.column_stack((xo[:, 0], yo[:, 0]))
    right = np.column_stack((xo[:, 1], yo[:, 1]))
    ring = bankline_polygon((left, right))
    return ring[:, 0], ring[:, 1]


def _courses(result, n_courses: int
             ) -> list[tuple[np.ndarray, np.ndarray]]:
    """Historical courses, evenly subsampled to at most *n_courses*.

    The final course is always kept: it is the one the eye reads first, and
    dropping it would leave the print without a subject.
    """
    history = list(getattr(result, "centerline_history", []) or [])
    line = result.centerline
    if not history:
        return [(line.x, line.y)]
    if len(history) > n_courses >= 1:
        index = np.unique(
            np.linspace(0, len(history) - 1, n_courses).round().astype(int))
        history = [history[k] for k in index]
    return history


def auto_window(result, widths: float = 45.0) -> tuple[float, float] | None:
    """A centred window of about *widths* channel widths of the reach.

    ``None`` if the reach is already that short, in which case there is nothing
    to crop. Used by the two bed styles, where the subject is the bar and pool
    of a bend rather than the shape of the belt.
    """
    line = result.centerline
    if line.length <= 0:
        return None
    span = widths * result.config.channel.width / line.length
    if span >= 1.0:
        return None
    return (0.5 - span / 2, 0.5 + span / 2)


def _window_mask(x: np.ndarray, y: np.ndarray, s: np.ndarray,
                 window: tuple[float, float] | None) -> np.ndarray:
    """Boolean mask of the nodes inside *window*, given as arc-length fractions."""
    if window is None or s.size == 0 or s[-1] <= 0:
        return np.ones(x.size, dtype=bool)
    start, end = window
    fraction = s / s[-1]
    return (fraction >= start) & (fraction <= end)


def _extent(result, style: ArtStyle) -> tuple[float, float, float, float]:
    """Bounding box of everything that will be drawn, plus the margin.

    With a window set, only the windowed stretch of the current centerline
    frames the print; the historical courses and oxbow lakes wander well beyond
    it and would pull the view back out to the whole belt.
    """
    line = result.centerline
    if style.window is not None:
        inside = _window_mask(line.x, line.y, line.s, style.window)
        if np.any(inside):
            x, y = line.x[inside], line.y[inside]
        else:
            x, y = line.x, line.y
    else:
        xs, ys = [line.x], [line.y]
        for cx, cy in _courses(result, style.n_courses):
            xs.append(np.asarray(cx))
            ys.append(np.asarray(cy))
        for cut in getattr(result, "cutoffs", []):
            if cut.oxbow_x.size:
                xs.append(cut.oxbow_x)
                ys.append(cut.oxbow_y)
        x, y = np.concatenate(xs), np.concatenate(ys)

    pad = style.margin * max(np.ptp(x), np.ptp(y), 1.0)
    half = result.config.channel.half_width
    return (x.min() - pad - half, x.max() + pad + half,
            y.min() - pad - half, y.max() + pad + half)


def _canvas(result, style: ArtStyle, ground: str):
    """A figure and one full-bleed axes with the reach framed and no decoration.

    The limits are widened here, once, to whichever of the two directions the
    canvas has spare room in, and the aspect is then locked with
    ``adjustable='box'``. Leaving it to ``adjustable='datalim'`` would change
    the limits *after* this function returned, so anything that read them --
    the blueprint grid, the scale bar, the vignette -- would be working from
    stale numbers and would cover only part of the canvas.
    """
    from matplotlib.figure import Figure

    x0, x1, y0, y1 = _extent(result, style)
    data_aspect = (x1 - x0) / max(y1 - y0, 1e-9)
    size = style.figsize(data_aspect)
    canvas_aspect = size[0] / size[1]

    # Grow the direction with room to spare, keeping the reach centred, so that
    # one data unit is the same length on both axes.
    if canvas_aspect > data_aspect:
        extra = (y1 - y0) * canvas_aspect - (x1 - x0)
        x0, x1 = x0 - extra / 2, x1 + extra / 2
    else:
        extra = (x1 - x0) / canvas_aspect - (y1 - y0)
        y0, y1 = y0 - extra / 2, y1 + extra / 2

    figure = Figure(figsize=size, dpi=style.dpi, facecolor=ground)
    axes = figure.add_axes((0.0, 0.0, 1.0, 1.0))
    axes.set_facecolor(ground)
    axes.set_xlim(x0, x1)
    axes.set_ylim(y0, y1)
    axes.set_aspect("equal", adjustable="box")
    axes.set_axis_off()
    return figure, axes


def _paper_texture(axes, style: ArtStyle, extent, strength: float,
                   colour: str = "#8A7A5C", zorder: int = -5) -> None:
    """Overlay a soft mottling so a flat fill reads as ink on paper.

    Two octaves of smoothed value noise, built with nothing but NumPy: a
    coarse blotch for the sheet and a finer grain for the tooth of the paper.
    """
    if strength <= 0:
        return
    from matplotlib.colors import LinearSegmentedColormap, to_rgb

    rng = np.random.default_rng(style.seed)
    field = np.zeros((256, 256))
    for octave, weight in ((8, 0.65), (32, 0.35)):
        coarse = rng.random((octave, octave))
        rows = np.linspace(0, octave - 1, 256)
        # Separable linear interpolation, cheap and smooth enough for texture.
        grid = np.interp(rows, np.arange(octave), np.arange(octave))
        lo = np.floor(grid).astype(int)
        hi = np.minimum(lo + 1, octave - 1)
        frac = (grid - lo)[:, None]
        vertical = coarse[lo] * (1 - frac) + coarse[hi] * frac
        frac = (grid - lo)[None, :]
        field += weight * (vertical[:, lo] * (1 - frac)
                           + vertical[:, hi] * frac)
    field -= field.min()
    field /= max(field.max(), 1e-9)

    red, green, blue = to_rgb(colour)
    cmap = LinearSegmentedColormap.from_list(
        "mottle", [(red, green, blue, 0.0), (red, green, blue, 1.0)])
    axes.imshow(field, extent=extent, origin="lower", cmap=cmap,
                alpha=0.10 * strength, interpolation="bilinear",
                zorder=zorder, aspect="auto")


def _vignette(axes, extent, colour: str = "#000000", strength: float = 0.25,
              zorder: int = 50) -> None:
    """Darken the corners, the way a printed sheet falls off at its edges."""
    if strength <= 0:
        return
    from matplotlib.colors import LinearSegmentedColormap, to_rgb

    grid = np.linspace(-1.0, 1.0, 256)
    xx, yy = np.meshgrid(grid, grid)
    radius = np.clip(np.hypot(xx, yy) / np.sqrt(2.0), 0.0, 1.0) ** 2.2
    red, green, blue = to_rgb(colour)
    cmap = LinearSegmentedColormap.from_list(
        "vig", [(red, green, blue, 0.0), (red, green, blue, 1.0)])
    axes.imshow(radius, extent=extent, origin="lower", cmap=cmap,
                alpha=strength, interpolation="bilinear", zorder=zorder,
                aspect="auto")


# --------------------------------------------------------------------------
# annotation
# --------------------------------------------------------------------------

def _default_title(result) -> str:
    """A title taken from the run: the river's name, or the curve's."""
    config = result.config
    name = (config.name or "").strip()
    if not name:
        if config.mode == "kinoshita":
            name = "Kinoshita Curve"
        else:
            name = Path(config.centerline_file).stem.replace("_", " ").title()
    return name


def _default_subtitle(result) -> str:
    """A subtitle stating what the print actually shows."""
    config = result.config
    steps = result.steps_completed
    if not steps:
        return "Synthetic riverbed — equilibrium bed topography"
    years = steps * config.migration.dt / 31_557_600.0
    span = (f"{years:,.0f} years" if years >= 2
            else f"{steps * config.migration.dt / 86400.0:,.0f} days")
    cutoffs = len(result.cutoffs)
    tail = f", {cutoffs} cutoff{'s' if cutoffs != 1 else ''}" if cutoffs else ""
    return f"Ancient courses over {span}{tail}"


def _cartouche(axes, style: ArtStyle, result, ink: str, ground: str,
               boxed: bool = True, loc: str = "lower left") -> None:
    """Title block, in the corner, restrained.

    Fisk's own title blocks are ornate; this is the same idea reduced to a
    rule and two lines of letterspaced type, which is what reads well at
    poster size without competing with the river.
    """
    title = style.title if style.title is not None else _default_title(result)
    if not title:
        return
    subtitle = (style.subtitle if style.subtitle is not None
                else _default_subtitle(result))

    corner = {"lower left": (0.045, 0.055, "left", "bottom"),
              "lower right": (0.955, 0.055, "right", "bottom"),
              "upper left": (0.045, 0.945, "left", "top"),
              "upper right": (0.955, 0.945, "right", "top")}[loc]
    px, py, ha, va = corner
    sign = -1.0 if va == "top" else 1.0

    box = dict(boxstyle="square,pad=0.9", facecolor=ground, edgecolor=ink,
               linewidth=0.8, alpha=0.88) if boxed else None
    text = axes.text(px, py + sign * 0.030, " ".join(title.upper()),
                     transform=axes.transAxes, ha=ha, va=va, color=ink,
                     fontsize=13, fontweight="bold", zorder=60)
    if box:
        text.set_bbox(box)
    if subtitle:
        axes.text(px, py, subtitle, transform=axes.transAxes, ha=ha, va=va,
                  color=ink, fontsize=8.5, alpha=0.85, style="italic",
                  zorder=60)


def _scale_bar(axes, style: ArtStyle, result, ink: str,
               loc: str = "lower right") -> None:
    """A scale bar in metres, rounded to a 1/2/5 x 10^n length.

    Sits in the opposite bottom corner from the cartouche by default, so the
    two never collide however wide the reach turns out to be.
    """
    if not style.scale_bar:
        return
    x0, x1 = axes.get_xlim()
    y0, y1 = axes.get_ylim()
    target = 0.15 * (x1 - x0)
    power = 10.0 ** np.floor(np.log10(max(target, 1.0)))
    length = next((step * power for step in (1, 2, 5, 10)
                   if step * power >= target * 0.6), power * 10)
    label = (f"{length / 1000:g} km" if length >= 1000 else f"{length:g} m")

    if loc == "lower right":
        bx = x1 - 0.045 * (x1 - x0) - length
    else:
        bx = x0 + 0.045 * (x1 - x0)
    by = y0 + 0.045 * (y1 - y0)
    thickness = 0.0055 * (y1 - y0)
    # Alternating filled and open cells, the surveyor's convention.
    cells = 4
    for cell in range(cells):
        axes.add_patch(_rect(bx + cell * length / cells, by,
                             length / cells, thickness,
                             ink if cell % 2 == 0 else "none", ink))
    axes.text(bx + length / 2, by + 2.2 * thickness, label, ha="center",
              va="bottom", color=ink, fontsize=8, zorder=60)


def _darken(colour: str, factor: float) -> tuple[float, float, float]:
    """Scale a colour towards black, for hatch lines and outlines."""
    from matplotlib.colors import to_rgb
    return tuple(channel * factor for channel in to_rgb(colour))


def _rect(x, y, width, height, face, edge):
    """A rectangle patch, kept local so the import stays lazy."""
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), width, height, facecolor=face, edgecolor=edge,
                     linewidth=0.6, zorder=60)


def _flow_arrow(axes, result, ink: str, alpha: float = 0.7) -> None:
    """A small arrow at the inlet, pointing the way the river flows."""
    line = result.centerline
    if line.n_nodes < 40:
        return
    step = max(line.n_nodes // 40, 1)
    x0, y0 = line.x[0], line.y[0]
    x1, y1 = line.x[step], line.y[step]
    scale = 2.5 * result.config.channel.width / max(np.hypot(x1 - x0, y1 - y0),
                                                    1e-9)
    axes.annotate("", xy=(x0 + (x1 - x0) * scale, y0 + (y1 - y0) * scale),
                  xytext=(x0, y0), zorder=60,
                  arrowprops=dict(arrowstyle="-|>", color=ink, alpha=alpha,
                                  linewidth=1.4, shrinkA=0, shrinkB=0))


# --------------------------------------------------------------------------
# styles
# --------------------------------------------------------------------------

def fisk(result, style: ArtStyle | None = None):
    """Harold Fisk's meander belt, drawn from a pyRiverBed run.

    Every historical course gets its own flat colour, oldest at the back, and
    the oxbow lakes are filled in the colour of the era that abandoned them.
    Older courses carry a hatch as well as a colour, which is how Fisk kept
    overlapping courses apart where their inks were close, and which keeps the
    print readable in greyscale.

    Parameters
    ----------
    result
        A finished :class:`~pyriverbed.model.RunResult`.
    style
        Print settings.

    Returns
    -------
    matplotlib.figure.Figure
    """
    style = style or ArtStyle()
    figure, axes = _canvas(result, style, FISK_PAPER)
    extent = (*axes.get_xlim(), *axes.get_ylim())
    half = result.config.channel.half_width

    _paper_texture(axes, style, extent, style.texture)

    courses = _courses(result, style.n_courses)
    n = len(courses)
    for age, (cx, cy) in enumerate(courses):
        newest = age == n - 1
        if newest and not style.show_current:
            continue
        colour = FISK_COURSES[age % len(FISK_COURSES)]
        hatch = "" if newest else FISK_HATCHES[age % len(FISK_HATCHES)]
        px, py = channel_ribbon(cx, cy, half)
        # The current course is Fisk's "mighty blank": paper, not ink, so the
        # live river reads as the hole punched through every older course.
        # Matplotlib draws hatching in the *edge* colour, so a hatched course
        # needs a darkened edge or the pattern is invisible against its own fill.
        axes.fill(px, py, facecolor=FISK_PAPER if newest else colour,
                  edgecolor=(FISK_INK if newest else
                             _darken(colour, 0.55) if hatch else colour),
                  linewidth=1.1 if newest else (0.5 if hatch else 0.0),
                  hatch=hatch or None, zorder=10 + age)

    if style.show_oxbows:
        cutoffs = sorted(getattr(result, "cutoffs", []), key=lambda c: c.step)
        total = max(result.steps_completed, 1)
        for cut in cutoffs:
            if not cut.oxbow_x.size:
                continue
            # Colour by era, matching the course that was live at the time.
            era = min(int(cut.step / total * max(n - 1, 1)), n - 1)
            colour = FISK_COURSES[era % len(FISK_COURSES)]
            px, py = channel_ribbon(cut.oxbow_x, cut.oxbow_y, half)
            axes.fill(px, py, facecolor=colour, edgecolor=colour,
                      linewidth=0.5, zorder=9)

    _flow_arrow(axes, result, FISK_INK)
    _scale_bar(axes, style, result, FISK_INK)
    _cartouche(axes, style, result, FISK_INK, FISK_PAPER)
    _vignette(axes, extent, "#5A4A32", 0.18)
    return figure


def strata(result, style: ArtStyle | None = None):
    """The channel courses as translucent layers in one perceptual ramp.

    Where :func:`fisk` separates the eras by hue, this separates them by
    lightness alone, so the print reads as a single object with depth in
    time rather than as a set of overlaid maps. Partial transparency lets the
    density of overprinting stand for how long the river stayed put.
    """
    style = style or ArtStyle()
    ground = "#141A1F"
    figure, axes = _canvas(result, style, ground)
    extent = (*axes.get_xlim(), *axes.get_ylim())
    half = result.config.channel.half_width

    import matplotlib.cm as cm
    from matplotlib.colors import Normalize

    courses = _courses(result, max(style.n_courses, 12))
    n = len(courses)
    ramp = cm.ScalarMappable(norm=Normalize(0, max(n - 1, 1)), cmap="magma")

    if style.show_oxbows:
        total = max(result.steps_completed, 1)
        for cut in sorted(getattr(result, "cutoffs", []), key=lambda c: c.step):
            if not cut.oxbow_x.size:
                continue
            px, py = channel_ribbon(cut.oxbow_x, cut.oxbow_y, half)
            axes.fill(px, py, color=ramp.to_rgba(cut.step / total * (n - 1)),
                      alpha=0.55, zorder=5)

    for age, (cx, cy) in enumerate(courses):
        px, py = channel_ribbon(cx, cy, half)
        axes.fill(px, py, color=ramp.to_rgba(age),
                  alpha=0.30 if age < n - 1 else 1.0, zorder=10 + age)

    _scale_bar(axes, style, result, "#D8D2C4")
    _cartouche(axes, style, result, "#EDE7DA", ground, boxed=False)
    _vignette(axes, extent, "#000000", 0.35)
    return figure


def blueprint(result, style: ArtStyle | None = None):
    """Cyanotype: a pale river on Prussian blue, with drafting furniture.

    The grid is spaced in channel widths rather than in metres, which is the
    scale the model actually thinks in, and it doubles as the print's texture.
    """
    style = style or ArtStyle()
    ground = "#0E3A5C"
    ink = "#D6E8F2"
    figure, axes = _canvas(result, style, ground)
    x0, x1, y0, y1 = *axes.get_xlim(), *axes.get_ylim()
    width = result.config.channel.width
    half = width / 2.0

    # Grid every 10 widths, with a hairline every width. Drawn with axvline and
    # axhline so it spans the whole canvas whatever the limits end up being.
    for spacing, alpha, lw in ((width, 0.07, 0.4), (10 * width, 0.18, 0.6)):
        for gx in np.arange(np.ceil(x0 / spacing) * spacing, x1, spacing):
            axes.axvline(gx, color=ink, alpha=alpha, linewidth=lw, zorder=1)
        for gy in np.arange(np.ceil(y0 / spacing) * spacing, y1, spacing):
            axes.axhline(gy, color=ink, alpha=alpha, linewidth=lw, zorder=1)

    if style.show_oxbows:
        for cut in getattr(result, "cutoffs", []):
            if not cut.oxbow_x.size:
                continue
            px, py = channel_ribbon(cut.oxbow_x, cut.oxbow_y, half)
            axes.fill(px, py, facecolor=ink, alpha=0.16, edgecolor=ink,
                      linewidth=0.5, zorder=8)

    for age, (cx, cy) in enumerate(_courses(result, style.n_courses)[:-1]):
        axes.plot(cx, cy, color=ink, alpha=0.22, linewidth=0.7, zorder=9)

    line = result.centerline
    px, py = channel_ribbon(line.x, line.y, half)
    axes.fill(px, py, facecolor=ink, alpha=0.92, edgecolor=ink,
              linewidth=0.8, zorder=12)
    axes.plot(line.x, line.y, color=ground, linewidth=0.5, alpha=0.5,
              linestyle=(0, (6, 4)), zorder=13)

    _flow_arrow(axes, result, ink)
    _scale_bar(axes, style, result, ink, loc="lower left")
    _cartouche(axes, style, result, ink, ground, boxed=True, loc="lower right")
    axes.text(0.045, 0.955,
              f"W = {width:g} m    H = {result.config.channel.depth:g} m    "
              f"grid = 1 W / 10 W",
              transform=axes.transAxes, ha="left", va="top", color=ink,
              fontsize=8, alpha=0.75, family="monospace", zorder=60)
    return figure


def nocturne(result, style: ArtStyle | None = None):
    """Warm metallic channel on near-black, with a soft glow.

    The glow is built by stroking the course several times at increasing width
    and decreasing alpha -- the cheapest way to fake a bloom, and it needs no
    image filter.
    """
    style = style or ArtStyle()
    ground = "#0B0C10"
    gold = "#E8C36A"
    figure, axes = _canvas(result, style, ground)
    extent = (*axes.get_xlim(), *axes.get_ylim())
    half = result.config.channel.half_width
    line = result.centerline

    if style.show_oxbows:
        for cut in getattr(result, "cutoffs", []):
            if not cut.oxbow_x.size:
                continue
            px, py = channel_ribbon(cut.oxbow_x, cut.oxbow_y, half)
            axes.fill(px, py, facecolor="#6E5A2E", alpha=0.45, zorder=5)
            axes.plot(px, py, color=gold, alpha=0.30, linewidth=0.6, zorder=6)

    for cx, cy in _courses(result, style.n_courses)[:-1]:
        axes.plot(cx, cy, color=gold, alpha=0.13, linewidth=0.6, zorder=8)

    for spread, alpha in ((9.0, 0.05), (6.0, 0.08), (3.5, 0.13), (1.8, 0.2)):
        axes.plot(line.x, line.y, color=gold, alpha=alpha,
                  linewidth=spread * 2.2, solid_capstyle="round", zorder=10)
    px, py = channel_ribbon(line.x, line.y, half)
    axes.fill(px, py, facecolor=gold, edgecolor="#FFF3D0", linewidth=0.6,
              zorder=14)

    _scale_bar(axes, style, result, "#9C8A63")
    _cartouche(axes, style, result, "#E8DFC8", ground, boxed=False)
    _vignette(axes, extent, "#000000", 0.45)
    return figure


def minimal(result, style: ArtStyle | None = None):
    """One ribbon, one accent colour, a lot of air.

    The whole point is restraint, so the oxbow lakes are drawn as unfilled
    outlines and the cartouche loses its box.
    """
    style = style or ArtStyle()
    ground = "#F4F1EA"
    accent = "#C2472F"
    figure, axes = _canvas(result, style, ground)
    half = result.config.channel.half_width
    line = result.centerline

    if style.show_oxbows:
        for cut in getattr(result, "cutoffs", []):
            if not cut.oxbow_x.size:
                continue
            px, py = channel_ribbon(cut.oxbow_x, cut.oxbow_y, half)
            axes.plot(px, py, color=accent, alpha=0.30, linewidth=0.7,
                      zorder=5)

    px, py = channel_ribbon(line.x, line.y, half)
    axes.fill(px, py, facecolor=accent, edgecolor="none", zorder=12)

    _cartouche(axes, style, result, "#2B2B2B", ground, boxed=False)
    if style.scale_bar:
        _scale_bar(axes, style, result, "#2B2B2B")
    return figure


def _bed_grid(result, n: int = 700, limits=None):
    """The bed point cloud resampled onto a regular grid, for contouring.

    Nearest-neighbour onto the grid via a k-d tree, masked outside the channel,
    which is enough for contour art and avoids pulling in a triangulator.

    Parameters
    ----------
    result
        A finished run.
    n
        Nominal grid resolution; the two directions are shared out between it
        according to the aspect ratio of the region.
    limits
        ``(x0, x1, y0, y1)`` to grid, or ``None`` for the whole cloud. Gridding
        only what is on the canvas is what keeps a zoomed bed print cheap: the
        cloud of a long migration run holds millions of points, almost all of
        them outside the window.
    """
    from scipy.spatial import cKDTree

    cloud = result.cloud
    x, y, z = cloud[:, 0], cloud[:, 1], cloud[:, 2]
    if limits is not None:
        x0, x1, y0, y1 = limits
        keep = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        if np.any(keep):
            x, y, z = x[keep], y[keep], z[keep]
    else:
        x0, x1, y0, y1 = x.min(), x.max(), y.min(), y.max()

    aspect = max(x1 - x0, 1e-9) / max(y1 - y0, 1e-9)
    nx = int(np.clip(n * np.sqrt(aspect), 80, 2400))
    ny = int(np.clip(n / max(np.sqrt(aspect), 1e-9), 80, 2400))
    mesh_x, mesh_y = np.meshgrid(np.linspace(x0, x1, nx),
                                 np.linspace(y0, y1, ny))

    tree = cKDTree(np.column_stack((x, y)))
    spacing = result.config.channel.interval
    distance, index = tree.query(np.column_stack((mesh_x.ravel(),
                                                  mesh_y.ravel())))
    grid = z[index].reshape(mesh_y.shape)
    # Outside the channel there is no data, so mask anything further from a
    # cloud point than the cloud's own spacing.
    grid = np.ma.masked_where(distance.reshape(mesh_y.shape) > 1.5 * spacing,
                              grid)
    return mesh_x, mesh_y, grid


def _bed_style(result, style: ArtStyle | None) -> ArtStyle:
    """Settings for a bed print: windowed, and captioned about the bed.

    The bed styles are the only ones whose subject is smaller than the reach,
    so they pick a window and a subtitle of their own unless the caller has
    already said what they want.
    """
    from dataclasses import replace

    style = style or ArtStyle()
    changes: dict = {}
    if style.window is None:
        window = auto_window(result)
        if window is not None:
            changes["window"] = window
    if style.subtitle is None:
        changes["subtitle"] = ("Synthetic bed topography — pools against the "
                              "outer banks, point bars against the inner")
    return replace(style, **changes) if changes else style


def bathymetry(result, style: ArtStyle | None = None):
    """The bed itself, as filled depth bands with crisp banklines.

    This is the carved-depth-chart look, and it is the one print that shows
    what pyRiverBed actually computes: the pools against the outer banks and
    the point bars against the inner ones. It zooms to a stretch of the reach
    for that reason -- see :attr:`ArtStyle.window`.
    """
    style = _bed_style(result, style)
    ground = "#F2EDE3"
    figure, axes = _canvas(result, style, ground)
    depth = result.config.channel.depth

    mesh_x, mesh_y, grid = _bed_grid(result, style.bed_grid,
                                     (*axes.get_xlim(), *axes.get_ylim()))
    levels = np.linspace(-1.05 * depth, 1.05 * depth, 15)
    axes.contourf(mesh_x, mesh_y, grid, levels=levels, cmap="YlGnBu_r",
                  extend="both", zorder=10)
    axes.contour(mesh_x, mesh_y, grid, levels=levels, colors="#1E3040",
                 linewidths=0.25, alpha=0.5, zorder=11)

    line = result.centerline
    px, py = channel_ribbon(line.x, line.y, result.config.channel.half_width)
    axes.plot(px, py, color="#12222E", linewidth=1.0, zorder=14)

    _flow_arrow(axes, result, "#12222E")
    _scale_bar(axes, style, result, "#12222E")
    _cartouche(axes, style, result, "#12222E", ground)
    return figure


def contour(result, style: ArtStyle | None = None):
    """The bed as thin contour lines only. Topographic line art.

    Zooms to a stretch of the reach, as :func:`bathymetry` does.
    """
    style = _bed_style(result, style)
    ground = "#FBF7F0"
    figure, axes = _canvas(result, style, ground)
    depth = result.config.channel.depth

    mesh_x, mesh_y, grid = _bed_grid(result, style.bed_grid,
                                     (*axes.get_xlim(), *axes.get_ylim()))
    axes.contour(mesh_x, mesh_y, grid,
                 levels=np.linspace(-depth, depth, 26), colors="#22333B",
                 linewidths=0.4, zorder=10)
    line = result.centerline
    px, py = channel_ribbon(line.x, line.y, result.config.channel.half_width)
    axes.plot(px, py, color="#22333B", linewidth=1.2, zorder=12)

    if style.show_oxbows:
        for cut in getattr(result, "cutoffs", []):
            if not cut.oxbow_x.size:
                continue
            ox, oy = channel_ribbon(cut.oxbow_x, cut.oxbow_y,
                                    result.config.channel.half_width)
            axes.plot(ox, oy, color="#22333B", alpha=0.35, linewidth=0.5,
                      zorder=8)

    _scale_bar(axes, style, result, "#22333B")
    _cartouche(axes, style, result, "#22333B", ground, boxed=False)
    return figure


#: Every style, by name. ``render`` and the ``pyriverbed art`` subcommand both
#: dispatch through this, so adding a style here is all it takes to expose it.
STYLES: dict[str, Callable] = {
    "fisk": fisk,
    "strata": strata,
    "blueprint": blueprint,
    "nocturne": nocturne,
    "minimal": minimal,
    "bathymetry": bathymetry,
    "contour": contour,
}


# --------------------------------------------------------------------------
# entry points
# --------------------------------------------------------------------------

def render(result, style_name: str = "fisk", style: ArtStyle | None = None,
           paths: Sequence[str | Path] = ()):
    """Draw one art print and optionally save it.

    Parameters
    ----------
    result
        A finished :class:`~pyriverbed.model.RunResult`.
    style_name
        Key into :data:`STYLES`.
    style
        Print settings. Defaults to :class:`ArtStyle`.
    paths
        Files to write. The extension picks the format, so ``.png`` for a
        raster print and ``.pdf`` or ``.svg`` for a vector one.

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    KeyError
        If *style_name* is not a known style.
    """
    if style_name not in STYLES:
        raise KeyError(f"unknown art style {style_name!r}; "
                       f"choose from {', '.join(sorted(STYLES))}")
    figure = STYLES[style_name](result, style or ArtStyle())
    for path in paths:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=(style or ArtStyle()).dpi,
                       facecolor=figure.get_facecolor())
        log.info("   wrote %s", path.name)
    return figure


def save_gallery(result, directory: str | Path, prefix: str = "art",
                 styles: Sequence[str] | None = None,
                 style: ArtStyle | None = None,
                 formats: Sequence[str] = ("png",)) -> list[Path]:
    """Render every style into *directory* and return the files written.

    Parameters
    ----------
    result
        A finished :class:`~pyriverbed.model.RunResult`.
    directory
        Where to write. Created if missing.
    prefix
        File name stem; the style name is appended.
    styles
        Which styles to render. ``None`` renders all of them.
    style
        Print settings, shared by every style.
    formats
        Extensions to write for each style.

    Returns
    -------
    list of Path
    """
    directory = Path(directory)
    chosen = list(styles) if styles else list(STYLES)
    written: list[Path] = []
    for name in chosen:
        paths = [directory / f"{prefix}_{name}.{ext}" for ext in formats]
        figure = render(result, name, style, paths)
        figure.clf()
        written += paths
    return written
