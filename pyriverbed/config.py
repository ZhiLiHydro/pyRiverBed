"""Configuration objects and input file handling for pyRiverBed.

The v1.x steering file was a bare column of 39 numbers whose meaning was
entirely positional: line 24 was the time step and there was no way to tell
that from looking at the file. v2.0 replaces it with a sectioned, commented
key/value file, and every value lands in a validated dataclass field.

Three input formats are understood, auto-detected by :func:`load_config`:

``*.ini``, ``*.cfg``, ``*.conf``
    The v2 format. Human-readable, self-documenting, order-independent.

``steering.txt`` (or any file that does not start with a section header)
    The v1.x positional format, still read so that old projects keep working.

Python
    Build :class:`Config` directly. This is what the notebook frontend does.
"""

from __future__ import annotations

import configparser
import math
from dataclasses import dataclass, fields, is_dataclass
from pathlib import Path
from typing import Any, Mapping

from ._version import __version__
from .logging_utils import get_logger

log = get_logger(__name__)

#: Elevations, curvatures and slopes below this are treated as exactly zero.
ZERO = 1e-8


class ConfigError(ValueError):
    """Raised when an input file is malformed or a parameter is out of range."""


# --------------------------------------------------------------------------
# sections
# --------------------------------------------------------------------------

@dataclass
class KinoshitaConfig:
    """Kinoshita curve parameters. Used by ``mode = kinoshita`` only.

    Attributes
    ----------
    n_bends
        Number of arc wavelengths to generate.
    arc_wavelength
        Arc wavelength (m).
    max_angular_amplitude
        Maximum angular amplitude (deg).
    skewness
        Skewness coefficient. Makes bends streamwise-asymmetric.
    flatness
        Flatness coefficient. Flattens bend apexes.
    """

    n_bends: int = 3
    arc_wavelength: float = 10.0
    max_angular_amplitude: float = 110.0
    skewness: float = 0.03125
    flatness: float = 0.00520833

    def validate(self) -> None:
        if self.n_bends < 1:
            raise ConfigError("kinoshita.n_bends must be at least 1")
        if self.arc_wavelength <= 0:
            raise ConfigError("kinoshita.arc_wavelength must be positive")

    @property
    def theta0(self) -> float:
        """Maximum angular amplitude in radians."""
        return math.radians(self.max_angular_amplitude)


@dataclass
class ChannelConfig:
    """Channel geometry and grid resolution.

    Attributes
    ----------
    width
        Channel width (m). Constant everywhere and always.
    depth
        Reach-averaged flow depth (m).
    slope
        Longitudinal channel slope. Set to 0 for a horizontal datum.
    transverse_slope_corrector
        Empirical multiplier on the transverse bed slope, in (0, 1]. The knob
        for calibrating bar-pool relief against measured bathymetry.
    ds
        Streamwise node spacing for the Kinoshita curve (m).
    n_offsets
        Number of polyline offsets per side of the centerline. Sets the
        transverse resolution, and through resampling the streamwise one too.
    """

    width: float = 0.6
    depth: float = 0.15
    slope: float = 0.0
    transverse_slope_corrector: float = 1.0
    ds: float = 0.03
    n_offsets: int = 10

    def validate(self) -> None:
        if self.width <= 0:
            raise ConfigError("channel.width must be positive")
        if self.depth <= 0:
            raise ConfigError("channel.depth must be positive")
        if self.slope < 0:
            raise ConfigError("channel.slope must not be negative")
        if not 0 < self.transverse_slope_corrector <= 1:
            raise ConfigError(
                "channel.transverse_slope_corrector must be in (0, 1]")
        if self.ds <= 0:
            raise ConfigError("channel.ds must be positive")
        if self.n_offsets < 1:
            raise ConfigError("channel.n_offsets must be at least 1")

    @property
    def half_width(self) -> float:
        """Channel half-width (m)."""
        return self.width / 2.0

    @property
    def interval(self) -> float:
        """Transverse node spacing (m)."""
        return self.half_width / self.n_offsets

    @property
    def beta(self) -> float:
        """Half-width-to-depth ratio."""
        return self.half_width / self.depth

    @property
    def scour_factor(self) -> float:
        """Beck (1988) scour factor A, from the half-width-to-depth ratio."""
        beta = self.beta
        return 3.8 * (1.0 + beta / 6.96 * math.exp(-6.96 / beta))


@dataclass
class LagConfig:
    """Curvature phase lag.

    Attributes
    ----------
    enabled
        Whether to phase-lag the curvature signal at all.
    strength
        Length of the upstream averaging window, in channel widths.
    """

    enabled: bool = True
    strength: float = 4.0

    def validate(self) -> None:
        if self.strength <= 0:
            raise ConfigError("lag.strength must be positive")


@dataclass
class CurvatureConfig:
    """Curvature estimation and smoothing.

    Attributes
    ----------
    method
        ``'arctan2'`` (default, signed), ``'cosine'`` (signed) or
        ``'circumcircle'`` (unsigned, for comparison only).
    smoothing_level
        Number of Savitzky-Golay passes applied **once**, when the centerline
        is built. This is the knob for cleaning up an imported centerline.
        Levels above 38 escalate as ``1.1 ** level`` to reach the very heavy
        smoothing that noisy satellite-derived centerlines sometimes need.
    migration_smoothing_level
        Number of passes applied **after every migration time step**, to remove
        the small-scale roughness that pointwise bank displacement introduces.

        This must be kept small, but not too small. Smoothing diffuses
        meanders, so a large value reapplied thousands of times destroys bends
        faster than the migration model can grow them and the reach decays
        towards a straight line. A single pass, on the other hand, does not
        keep up with the node-scale roughness: about 1% of nodes end up
        carrying curvature spikes an order of magnitude beyond anything
        physical, which inflate the sinuosity with wiggles that are not bends.

        Two passes is the knee, and the outcome is insensitive between 2 and
        about 8, so 2 is the default.

        v1.x had a single smoothing level and used it in both places, so a
        migration run with the GUI default of 20 flattened itself. Splitting
        the two is what makes long runs behave.
    despike
        Whether to replace isolated curvature spikes by the neighbour mean.
    """

    method: str = "arctan2"
    smoothing_level: int = 20
    migration_smoothing_level: int = 2
    despike: bool = True

    METHODS = ("arctan2", "cosine", "circumcircle")

    def validate(self) -> None:
        if self.method not in self.METHODS:
            raise ConfigError(
                f"curvature.method must be one of {self.METHODS}, "
                f"got {self.method!r}")
        if self.smoothing_level < 0:
            raise ConfigError("curvature.smoothing_level must not be negative")
        if self.migration_smoothing_level < 0:
            raise ConfigError(
                "curvature.migration_smoothing_level must not be negative")

    @staticmethod
    def _passes(level: int) -> int:
        """Convert a smoothing level to a number of filter passes."""
        return level if level < 39 else int(round(1.1 ** level))

    @property
    def n_passes(self) -> int:
        """Smoothing passes applied when the centerline is built."""
        return self._passes(self.smoothing_level)

    @property
    def n_migration_passes(self) -> int:
        """Smoothing passes applied after each migration time step."""
        return self._passes(self.migration_smoothing_level)


@dataclass
class FlipConfig:
    """Reflections applied to the planform and to the bed.

    Attributes
    ----------
    streamwise
        Reverse the flow direction of the centerline.
    transverse
        Mirror the bed about the centerline. Purely cosmetic: the default of
        ``no`` already puts the pool on the outer bank. (In v1.x this had to be
        switched on to get the pool on the right side, because the transverse
        slope carried the wrong sign; see :func:`pyriverbed.bed.transverse_slope`.)
    """

    streamwise: bool = False
    transverse: bool = False

    def validate(self) -> None:
        return None


@dataclass
class OutputConfig:
    """What to write and where.

    Attributes
    ----------
    directory
        Output directory. Created if missing.
    save_xyz
        Write the riverbed point cloud (``*.xyz``).
    save_bankline
        Write the banklines as a closed polyline (``*.i2s``).
    save_mesh
        Write the finite element mesh and boundary conditions
        (``*.t3s``, ``*.dat``, ``*.cli``, ``*.bc2``).
    save_figures
        Write the summary figure of each graphic printout.
    figure_formats
        Formats for the figure written at time step 0.
    save_gif
        Assemble the per-step frames into animated GIFs. Needs ``imageio``.
    gif_fps
        Frames per second of the GIFs.
    frame_dpi
        Resolution of the per-step frames.
    figure_dpi
        Resolution of the headline figure.
    log_file
        Name of the run log inside ``directory``. Empty disables file logging.
    """

    directory: str = "."
    save_xyz: bool = True
    save_bankline: bool = True
    save_mesh: bool = True
    save_figures: bool = True
    figure_formats: tuple[str, ...] = ("png", "pdf")
    save_gif: bool = True
    gif_fps: int = 24
    frame_dpi: int = 150
    figure_dpi: int = 300
    log_file: str = "pyriverbed.log"

    def validate(self) -> None:
        if self.gif_fps < 1:
            raise ConfigError("output.gif_fps must be at least 1")
        if self.frame_dpi < 10 or self.figure_dpi < 10:
            raise ConfigError("output DPI values must be at least 10")

    @property
    def path(self) -> Path:
        """Output directory as a :class:`~pathlib.Path`."""
        return Path(self.directory)


@dataclass
class MigrationConfig:
    """Meander migration, after Ikeda, Parker & Sawai (1981).

    Attributes
    ----------
    enabled
        Whether to migrate the planform at all. With migration off, pyRiverBed
        generates one equilibrium bed for the given planform and stops.
    n_steps
        Number of time steps.
    dt
        Time step (s).
    e0
        Bank erosion coefficient (1/s). The single rate-setting parameter.
    ub0
        Amplitude of the random inlet velocity perturbation. Must be non-zero
        to grow meanders out of an initially straight channel.
    c0
        Deterministic inlet velocity perturbation.
    cf0
        Reach-averaged friction coefficient.
    fr0
        Reach-averaged Froude number.
    end_taper_widths
        Length, in channel widths, over which the bank displacement is ramped
        from zero at each end of the reach to its full value in the interior.

        The linear theory has no upstream reach to spread an inlet perturbation
        over, so with a free inlet the ``ub0`` noise and the curvature feedback
        amplify each other into a sharp hook within the first few widths --
        curvature far outside the range where the theory is valid, which then
        coils up and fires spurious neck cutoffs. Holding the end nodes fixed
        is the standard remedy and leaves the interior physics untouched.

        The default of 2 widths covers the inlet transient, whose length is
        ``1 / (2 cf0 beta chi)``. Set to 0 for the free ends of v1.x, which are
        only stable under heavy per-step smoothing.
    log_every
        Time step interval for progress logging.
    plot_every
        Time step interval for figures and animation frames.
    seed
        Seed for the random number generator. ``None`` leaves it unseeded;
        set an integer to make a stochastic run reproducible.
    """

    enabled: bool = False
    n_steps: int = 10000
    dt: float = 86400.0
    e0: float = 1e-7
    ub0: float = 0.0
    c0: float = 0.0
    cf0: float = 0.01
    fr0: float = 0.1
    end_taper_widths: float = 2.0
    log_every: int = 50
    plot_every: int = 100
    seed: int | None = None

    def validate(self) -> None:
        if self.n_steps < 0:
            raise ConfigError("migration.n_steps must not be negative")
        if self.dt <= 0:
            raise ConfigError("migration.dt must be positive")
        if self.cf0 <= 0:
            raise ConfigError("migration.cf0 must be positive")
        if self.fr0 <= 0:
            raise ConfigError("migration.fr0 must be positive")
        if self.end_taper_widths < 0:
            raise ConfigError(
                "migration.end_taper_widths must not be negative")
        if self.log_every < 1:
            raise ConfigError("migration.log_every must be at least 1")
        if self.plot_every < 1:
            raise ConfigError("migration.plot_every must be at least 1")


@dataclass
class NeckCutoffConfig:
    """Neck cutoff detection.

    A neck cutoff is declared for the first pair of nodes that are far apart
    along the channel but close together in space.

    Attributes
    ----------
    enabled
        Whether to detect neck cutoffs.
    min_separation_widths
        Minimum along-channel separation of the two nodes, in channel widths.
    max_distance_widths
        Maximum straight-line distance between the two nodes, in channel
        widths.
    end_margin_widths
        Length, in channel widths, kept clear of each end of the reach. The
        ends carry the artificial straight extensions and the inlet transient
        of the migration model, so a cutoff detected there is a boundary
        artifact rather than a meander closing on itself. Set to 0 to search
        the whole reach, as v1.x did.
    """

    enabled: bool = True
    min_separation_widths: float = 2.0
    max_distance_widths: float = 1.0
    end_margin_widths: float = 2.0

    def validate(self) -> None:
        if self.min_separation_widths <= 0:
            raise ConfigError(
                "neck_cutoff.min_separation_widths must be positive")
        if self.max_distance_widths <= 0:
            raise ConfigError(
                "neck_cutoff.max_distance_widths must be positive")
        if self.end_margin_widths < 0:
            raise ConfigError(
                "neck_cutoff.end_margin_widths must not be negative")


@dataclass
class ChuteCutoffConfig:
    """Chute cutoff modeling.

    Chute cutoffs are *conditionally random*: the planform geometry decides
    where a chute channel is possible, and a random draw decides whether one
    is actually carved. See ``THEORY_GUIDE.md`` for the reasoning.

    Attributes
    ----------
    enabled
        Whether to model chute cutoffs. Needs migration to be enabled.
    frequency
        Probability of triggering a chute cutoff in one time step, in [0, 1].
        ``0.1`` means a 10% chance per step, i.e. a mean recurrence of ten
        steps. Calibrate as ``dt / recurrence_interval``.
    start_step
        Time steps of spin-up before chute cutoffs are allowed, so that bends
        have developed before any of them can be bypassed.
    entrance
        ``'apex'`` to start and end chute channels at bend apexes,
        ``'inflection'`` to use the curvature inflection points.
    span
        Number of entrance points a chute channel spans. ``2`` bypasses one
        full meander loop.
    max_valley_angle
        Maximum angle between the chute chord and the valley axis (deg). A
        small angle means a large slope advantage. ``90`` disables the test.
    min_length_widths
        Minimum length of the bypassed channel reach, in channel widths.
    min_sinuosity
        Minimum sinuosity of the bypassed reach, i.e. its arc length divided by
        the chute chord length. This is the slope advantage the chute would
        gain, so it is the most direct statement of why a chute forms at all.
        Values at or below 1 admit chutes across essentially straight reaches,
        which is unphysical.
    end_margin
        Entrance points kept away from each end of the centerline, where the
        planform is contaminated by the straight extensions and by the inlet
        transient.
    """

    enabled: bool = False
    frequency: float = 0.1
    start_step: int = 2000
    entrance: str = "apex"
    span: int = 2
    max_valley_angle: float = 30.0
    min_length_widths: float = 10.0
    min_sinuosity: float = 1.2
    end_margin: int = 3

    ENTRANCES = ("apex", "inflection")

    def validate(self) -> None:
        if not 0.0 <= self.frequency <= 1.0:
            raise ConfigError(
                "chute_cutoff.frequency must be a probability in [0, 1]")
        if self.start_step < 0:
            raise ConfigError("chute_cutoff.start_step must not be negative")
        if self.entrance not in self.ENTRANCES:
            raise ConfigError(
                f"chute_cutoff.entrance must be one of {self.ENTRANCES}, "
                f"got {self.entrance!r}")
        if self.span < 1:
            raise ConfigError("chute_cutoff.span must be at least 1")
        if not 0.0 <= self.max_valley_angle <= 90.0:
            raise ConfigError(
                "chute_cutoff.max_valley_angle must be in [0, 90] degrees")
        if self.min_length_widths <= 0:
            raise ConfigError(
                "chute_cutoff.min_length_widths must be positive")
        if self.min_sinuosity < 1.0:
            raise ConfigError(
                "chute_cutoff.min_sinuosity must be at least 1 (a reach "
                "cannot be shorter than the straight line across it)")
        if self.end_margin < 0:
            raise ConfigError("chute_cutoff.end_margin must not be negative")


# --------------------------------------------------------------------------
# top level
# --------------------------------------------------------------------------

@dataclass
class Config:
    """Complete pyRiverBed configuration.

    Attributes
    ----------
    mode
        ``'kinoshita'`` to synthesise a Kinoshita curve, ``'centerline'`` to
        read a centerline from ``centerline_file``.
    centerline_file
        Two-column text file of centerline coordinates, used by
        ``mode = centerline``.
    name
        Stem used for every output file. Empty derives it from the mode or
        from the centerline file name.
    """

    mode: str = "kinoshita"
    centerline_file: str = "mycenterline.txt"
    name: str = ""
    kinoshita: KinoshitaConfig = None  # type: ignore[assignment]
    channel: ChannelConfig = None  # type: ignore[assignment]
    curvature: CurvatureConfig = None  # type: ignore[assignment]
    lag: LagConfig = None  # type: ignore[assignment]
    flip: FlipConfig = None  # type: ignore[assignment]
    output: OutputConfig = None  # type: ignore[assignment]
    migration: MigrationConfig = None  # type: ignore[assignment]
    neck_cutoff: NeckCutoffConfig = None  # type: ignore[assignment]
    chute_cutoff: ChuteCutoffConfig = None  # type: ignore[assignment]

    MODES = ("kinoshita", "centerline")

    #: Section name -> dataclass, in the order they appear in an input file.
    SECTIONS: Mapping[str, type] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        for name, factory in (
            ("kinoshita", KinoshitaConfig),
            ("channel", ChannelConfig),
            ("curvature", CurvatureConfig),
            ("lag", LagConfig),
            ("flip", FlipConfig),
            ("output", OutputConfig),
            ("migration", MigrationConfig),
            ("neck_cutoff", NeckCutoffConfig),
            ("chute_cutoff", ChuteCutoffConfig),
        ):
            if getattr(self, name) is None:
                setattr(self, name, factory())

    # -- derived ----------------------------------------------------------

    @property
    def stem(self) -> str:
        """Stem of every output file name."""
        if self.name:
            return self.name
        if self.mode == "kinoshita":
            return "kinoshita"
        return Path(self.centerline_file).stem or "centerline"

    @property
    def n_steps(self) -> int:
        """Number of migration time steps actually run."""
        return self.migration.n_steps if self.migration.enabled else 0

    @property
    def chute_active(self) -> bool:
        """Whether chute cutoffs can happen in this run."""
        return self.chute_cutoff.enabled and self.migration.enabled

    @property
    def neck_active(self) -> bool:
        """Whether neck cutoffs can happen in this run."""
        return self.neck_cutoff.enabled and self.migration.enabled

    # -- validation -------------------------------------------------------

    def validate(self) -> "Config":
        """Validate every section. Returns ``self`` so calls can be chained.

        Raises
        ------
        ConfigError
            If any parameter is out of range or inconsistent.
        """
        if self.mode not in self.MODES:
            raise ConfigError(
                f"mode must be one of {self.MODES}, got {self.mode!r}")
        for f in fields(self):
            value = getattr(self, f.name)
            if is_dataclass(value):
                value.validate()
        if self.mode == "centerline" and not self.centerline_file:
            raise ConfigError(
                "mode = centerline needs centerline_file to be set")
        if self.chute_cutoff.enabled and not self.migration.enabled:
            log.warning(
                "chute cutoffs are enabled but migration is off; "
                "no cutoff can happen in a run with no time stepping")
        if self.migration.enabled:
            step = self.migration.e0 * self.migration.dt
            if step > 0.1:
                log.warning(
                    "e0 * dt = %.3g channel widths per step is large; the "
                    "planform may become unstable. Consider reducing e0 or dt",
                    step)
        return self

    # -- conversion -------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """Return a nested plain-``dict`` view of the configuration."""
        out: dict[str, Any] = {}
        for f in fields(self):
            value = getattr(self, f.name)
            if is_dataclass(value):
                out[f.name] = {sf.name: getattr(value, sf.name)
                               for sf in fields(value)}
            else:
                out[f.name] = value
        return out

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Config":
        """Build a :class:`Config` from a nested mapping.

        Unknown keys raise :class:`ConfigError` rather than being ignored, so
        a typo in an input file is reported instead of silently doing nothing.
        """
        section_types = {
            "kinoshita": KinoshitaConfig,
            "channel": ChannelConfig,
            "curvature": CurvatureConfig,
            "lag": LagConfig,
            "flip": FlipConfig,
            "output": OutputConfig,
            "migration": MigrationConfig,
            "neck_cutoff": NeckCutoffConfig,
            "chute_cutoff": ChuteCutoffConfig,
        }
        scalar_names = {"mode", "centerline_file", "name"}
        kwargs: dict[str, Any] = {}
        for key, value in data.items():
            if key in scalar_names:
                kwargs[key] = value
            elif key in section_types:
                kwargs[key] = _build_section(section_types[key], value, key)
            else:
                raise ConfigError(f"unknown configuration key {key!r}")
        return cls(**kwargs)


def _build_section(section_type: type, values: Mapping[str, Any],
                   section_name: str):
    """Instantiate one configuration section, coercing every value."""
    if not isinstance(values, Mapping):
        raise ConfigError(f"section [{section_name}] must be a mapping")
    known = {f.name: f for f in fields(section_type)}
    kwargs = {}
    for key, raw in values.items():
        if key not in known:
            raise ConfigError(
                f"unknown key {key!r} in section [{section_name}]; "
                f"valid keys are {', '.join(sorted(known))}")
        kwargs[key] = _coerce(raw, known[key], section_name, key)
    return section_type(**kwargs)


_TRUE = {"1", "true", "yes", "on", "y", "t"}
_FALSE = {"0", "false", "no", "off", "n", "f"}


def _coerce(raw: Any, dc_field, section: str, key: str) -> Any:
    """Coerce *raw* to the declared type of a dataclass field."""
    target = dc_field.type
    if isinstance(target, str):  # postponed annotations
        target = target.replace("int | None", "optional_int")
    where = f"[{section}] {key}"

    def as_bool(value: Any) -> bool:
        if isinstance(value, bool):
            return value
        text = str(value).strip().lower()
        if text in _TRUE:
            return True
        if text in _FALSE:
            return False
        raise ConfigError(f"{where}: expected a yes/no value, got {value!r}")

    try:
        if target in (bool, "bool"):
            return as_bool(raw)
        if target in (int, "int"):
            return int(float(str(raw).strip()))
        if target in (float, "float"):
            return float(str(raw).strip())
        if target in (str, "str"):
            return str(raw).strip()
        if target == "optional_int":
            text = str(raw).strip().lower()
            if isinstance(raw, int):
                return raw
            if text in ("", "none", "null", "-"):
                return None
            return int(float(text))
        if target in ("tuple[str, ...]",):
            if isinstance(raw, (list, tuple)):
                return tuple(str(v).strip() for v in raw)
            return tuple(v.strip() for v in str(raw).replace(",", " ").split()
                         if v.strip())
    except ConfigError:
        raise
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{where}: cannot read {raw!r} as {target} "
                          f"({exc})") from exc
    return raw


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def load_config(path: str | Path) -> Config:
    """Read a configuration from *path*, detecting the format.

    Parameters
    ----------
    path
        An ``*.ini``-style v2 input file, or a v1.x positional
        ``steering.txt``.

    Returns
    -------
    Config
        The validated configuration.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ConfigError
        If the file cannot be parsed.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"input file not found: {path}")
    text = path.read_text(encoding="utf-8")
    if _looks_like_ini(text):
        log.debug("reading %s as a v2 input file", path)
        return read_ini(text)
    log.warning("%s looks like a v1.x steering file; reading it in legacy "
                "mode. Run 'pyriverbed convert' to migrate it to the v2 "
                "format", path)
    return read_legacy_steering(text)


def _looks_like_ini(text: str) -> bool:
    """Whether *text* starts with an INI section header."""
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped[0] in "#;":
            continue
        return stripped.startswith("[")
    return False


def read_ini(text: str) -> Config:
    """Parse the v2 input format from a string."""
    parser = configparser.ConfigParser(
        inline_comment_prefixes=("#", ";"),
        interpolation=None,
    )
    parser.optionxform = str  # keep key case as written
    try:
        parser.read_string(text)
    except configparser.Error as exc:
        raise ConfigError(f"cannot parse input file: {exc}") from exc

    data: dict[str, Any] = {}
    for section in parser.sections():
        key = section.strip().lower().replace(" ", "_").replace("-", "_")
        items = {k.strip().lower(): v for k, v in parser.items(section)}
        if key in ("run", "general", "pyriverbed"):
            # Top-level scalars live in a section so the file stays valid INI.
            for scalar in ("mode", "centerline_file", "name"):
                if scalar in items:
                    data[scalar] = items[scalar].strip()
            unknown = set(items) - {"mode", "centerline_file", "name",
                                    "version"}
            if unknown:
                raise ConfigError(
                    f"unknown key(s) in section [{section}]: "
                    f"{', '.join(sorted(unknown))}")
        else:
            data[key] = items
    return Config.from_dict(data).validate()


#: Order of the 39 values in a v1.x steering file, as read by v1.1.0.
_LEGACY_ORDER: tuple[tuple[str, str], ...] = (
    ("mode", ""),
    ("kinoshita", "n_bends"),
    ("kinoshita", "arc_wavelength"),
    ("kinoshita", "max_angular_amplitude"),
    ("kinoshita", "skewness"),
    ("kinoshita", "flatness"),
    ("channel", "width"),
    ("channel", "depth"),
    ("channel", "slope"),
    ("channel", "ds"),
    ("channel", "n_offsets"),
    ("lag", "enabled"),
    ("lag", "strength"),
    ("output", "save_xyz"),
    ("output", "save_bankline"),
    ("output", "save_mesh"),
    ("flip", "streamwise"),
    ("flip", "transverse"),
    ("migration", "enabled"),
    ("migration", "ub0"),
    ("migration", "c0"),
    ("migration", "cf0"),
    ("migration", "fr0"),
    ("migration", "dt"),
    ("migration", "e0"),
    ("migration", "log_every"),
    ("migration", "n_steps"),
    ("migration", "plot_every"),
    ("output", "gif_fps"),
    ("curvature", "smoothing_level"),
    ("channel", "transverse_slope_corrector"),
    ("chute_cutoff", "enabled"),
    ("chute_cutoff", "frequency"),
    ("chute_cutoff", "start_step"),
    ("chute_cutoff", "entrance"),
    ("chute_cutoff", "span"),
    ("chute_cutoff", "max_valley_angle"),
    ("chute_cutoff", "min_length_widths"),
    ("chute_cutoff", "end_margin"),
)


def read_legacy_steering(text: str) -> Config:
    """Parse a v1.x ``steering.txt``.

    The first line is the centerline file name, followed by one number per
    line in the fixed v1.x order. Files written before v1.1.0 stop after 31
    numbers; the chute cutoff parameters then keep their defaults.
    """
    lines = [line.strip() for line in text.splitlines()]
    lines = [line for line in lines if line != ""]
    if not lines:
        raise ConfigError("steering file is empty")
    centerline_file = lines[0]
    values = lines[1:]
    if len(values) < 31:
        raise ConfigError(
            f"steering file has {len(values)} values, expected at least 31")

    data: dict[str, Any] = {"centerline_file": centerline_file}
    sections: dict[str, dict[str, Any]] = {}
    for raw, (section, key) in zip(values, _LEGACY_ORDER):
        if section == "mode":
            data["mode"] = "kinoshita" if int(float(raw)) == 1 else "centerline"
            continue
        if (section, key) == ("chute_cutoff", "entrance"):
            raw = "apex" if int(float(raw)) == 1 else "inflection"
        elif (section, key) == ("flip", "transverse"):
            # v1.x had the transverse bed slope sign inverted and cancelled it
            # with this flip. v2 fixes the sign, so the flag has to be inverted
            # for a v1.x file to reproduce v1.x output.
            raw = 0 if int(float(raw)) != 0 else 1
        sections.setdefault(section, {})[key] = raw
    data.update(sections)
    config = Config.from_dict(data)
    if len(values) < len(_LEGACY_ORDER):
        log.info("steering file predates v1.1.0; chute cutoff modeling is off")
        config.chute_cutoff = ChuteCutoffConfig()
    # v1.x used one smoothing level both for building the centerline and for
    # every migration step. Reproduce that here so a legacy file gives legacy
    # results, and warn if it is a value that will flatten the reach.
    config.curvature.migration_smoothing_level = config.curvature.smoothing_level
    # v1.x migrated its end nodes freely and searched the whole reach for neck
    # cutoffs. Both are only stable under v1.x's heavy per-step smoothing, but
    # switching them off here is what makes a legacy file reproduce legacy
    # output.
    config.migration.end_taper_widths = 0.0
    config.neck_cutoff.end_margin_widths = 0.0
    if config.migration.enabled and config.curvature.smoothing_level > 3:
        log.warning(
            "this v1.x file smooths %d times per migration step, which "
            "diffuses meanders faster than they grow. Set "
            "curvature.migration_smoothing_level = 2 in the converted v2 file",
            config.curvature.n_passes)
    return config.validate()


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------

#: One short comment per key, emitted above it when writing an input file.
_COMMENTS: Mapping[str, str] = {
    "mode": "kinoshita = synthesise a curve | centerline = read from file",
    "centerline_file": "two-column x y text file, used by mode = centerline",
    "name": "output file stem; empty = derive from mode / file name",
    "kinoshita.n_bends": "number of arc wavelengths to generate",
    "kinoshita.arc_wavelength": "arc wavelength (m)",
    "kinoshita.max_angular_amplitude": "maximum angular amplitude (deg)",
    "kinoshita.skewness": "streamwise asymmetry of the bends",
    "kinoshita.flatness": "flattening of the bend apexes",
    "channel.width": "channel width (m), constant everywhere",
    "channel.depth": "reach-averaged flow depth (m)",
    "channel.slope": "longitudinal channel slope; 0 = horizontal datum",
    "channel.transverse_slope_corrector":
        "(0, 1] multiplier on transverse bed slope; lower = milder bar-pool "
        "relief",
    "channel.ds": "streamwise node spacing (m), mode = kinoshita only",
    "channel.n_offsets": "polyline offsets per side; sets the resolution",
    "curvature.method": "arctan2 | cosine | circumcircle",
    "curvature.smoothing_level":
        "Savitzky-Golay passes applied once when building the centerline; the "
        "most important number for a real river",
    "curvature.migration_smoothing_level":
        "passes after each migration step; 2-8 works, below that node-scale "
        "spikes survive, above it meanders diffuse away",
    "curvature.despike": "replace isolated curvature spikes",
    "lag.enabled": "phase-lag the curvature signal",
    "lag.strength": "upstream averaging window, in channel widths",
    "flip.streamwise": "reverse the flow direction",
    "flip.transverse": "mirror the bed about the centerline",
    "output.directory": "output directory, created if missing",
    "output.save_xyz": "riverbed point cloud (.xyz)",
    "output.save_bankline": "banklines as a closed polyline (.i2s)",
    "output.save_mesh": "FEM mesh and BC files (.t3s .dat .cli .bc2)",
    "output.save_figures": "summary figure at every graphic printout",
    "output.figure_formats": "formats of the headline figure",
    "output.save_gif": "assemble frames into GIFs (needs imageio)",
    "output.gif_fps": "frames per second of the GIFs",
    "output.frame_dpi": "resolution of the animation frames",
    "output.figure_dpi": "resolution of the headline figure",
    "output.log_file": "run log inside the output directory; empty = none",
    "migration.enabled": "migrate the planform in time",
    "migration.n_steps": "number of time steps",
    "migration.dt": "time step (s)",
    "migration.e0": "bank erosion coefficient (1/s); the rate-setting knob",
    "migration.ub0": "random inlet perturbation; needed to grow meanders "
                     "from a straight channel",
    "migration.c0": "deterministic inlet perturbation",
    "migration.cf0": "reach-averaged friction coefficient",
    "migration.fr0": "reach-averaged Froude number",
    "migration.end_taper_widths":
        "ramp bank displacement to zero over this many widths at each end; "
        "keeps the free-inlet instability out. 0 = v1.x free ends",
    "migration.log_every": "time steps between progress reports",
    "migration.plot_every": "time steps between figures / frames",
    "migration.seed": "RNG seed; empty = unseeded",
    "neck_cutoff.enabled": "detect neck cutoffs",
    "neck_cutoff.min_separation_widths":
        "minimum along-channel separation of the two nodes (channel widths)",
    "neck_cutoff.max_distance_widths":
        "maximum straight-line distance between them (channel widths)",
    "neck_cutoff.end_margin_widths":
        "reach ends kept clear of the search (channel widths); cutoffs there "
        "are boundary artifacts. 0 = search everything, as v1.x did",
    "chute_cutoff.enabled": "model chute cutoffs (needs migration)",
    "chute_cutoff.frequency":
        "probability per time step; 0.1 = 10% chance each step",
    "chute_cutoff.start_step": "spin-up before cutoffs are allowed",
    "chute_cutoff.entrance": "apex | inflection",
    "chute_cutoff.span": "entrance points spanned; 2 = one meander loop",
    "chute_cutoff.max_valley_angle":
        "max angle with the valley axis (deg); 90 = no alignment test",
    "chute_cutoff.min_length_widths":
        "minimum bypassed reach length (channel widths)",
    "chute_cutoff.min_sinuosity":
        "minimum sinuosity of the bypassed reach = the chute's slope advantage",
    "chute_cutoff.end_margin": "entrance points kept clear of the reach ends",
}

_SECTION_HEADERS: Mapping[str, str] = {
    "run": "Which planform to work on",
    "channel": "Channel geometry and grid resolution",
    "kinoshita": "Kinoshita curve, used by mode = kinoshita only",
    "curvature": "Curvature estimation and centerline smoothing",
    "lag": "Curvature phase lag (secondary-flow memory)",
    "flip": "Reflections",
    "output": "What to write and where",
    "migration": "Meander migration (Ikeda, Parker & Sawai, 1981)",
    "neck_cutoff": "Neck cutoffs: geometric, deterministic",
    "chute_cutoff": "Chute cutoffs: geometric criteria + a random trigger",
}

_SECTION_ORDER = ("run", "channel", "kinoshita", "curvature", "lag", "flip",
                  "output", "migration", "neck_cutoff", "chute_cutoff")


def _format_value(value: Any) -> str:
    """Render a value the way the input file should show it."""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if value is None:
        return ""
    if isinstance(value, (tuple, list)):
        return ", ".join(str(v) for v in value)
    if isinstance(value, float):
        if value == 0:
            return "0"
        if 1e-4 <= abs(value) < 1e6:
            return f"{value:g}"
        return f"{value:.6g}"
    return str(value)


def dump_ini(config: Config, header: bool = True) -> str:
    """Render *config* as a commented v2 input file.

    Parameters
    ----------
    config
        The configuration to write.
    header
        Whether to include the explanatory file header.

    Returns
    -------
    str
        The file contents.
    """
    data = config.to_dict()
    lines: list[str] = []
    if header:
        lines += [
            f"# pyRiverBed {__version__} input file",
            "#",
            "# Sections and keys may appear in any order. Values are",
            "# yes/no for switches, plain numbers otherwise. Anything after",
            "# a '#' is a comment. See THEORY_GUIDE.md for the physics and",
            "# README.md for how to run.",
            "",
        ]
    for section in _SECTION_ORDER:
        lines.append(f"# --- {_SECTION_HEADERS[section]} "
                     + "-" * max(0, 66 - len(_SECTION_HEADERS[section])))
        lines.append(f"[{section}]")
        if section == "run":
            items = [(k, data[k]) for k in ("mode", "centerline_file", "name")]
            prefix = ""
        else:
            items = list(data[section].items())
            prefix = f"{section}."
        width = max(len(k) for k, _ in items)
        for key, value in items:
            comment = _COMMENTS.get(f"{prefix}{key}")
            rendered = _format_value(value)
            line = f"{key.ljust(width)} = {rendered}"
            if comment:
                line = f"{line.ljust(width + 26)}  # {comment}"
            lines.append(line)
        lines.append("")
    return "\n".join(lines).rstrip("\n") + "\n"


def write_config(config: Config, path: str | Path) -> Path:
    """Write *config* to *path* in the v2 format.

    Returns
    -------
    pathlib.Path
        The path written.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dump_ini(config), encoding="utf-8")
    log.info("wrote input file %s", path)
    return path


def default_config() -> Config:
    """Return the default configuration.

    The defaults reproduce the laboratory flume of Abad & Garcia (2009), which
    is the worked example in the README.
    """
    return Config().validate()
