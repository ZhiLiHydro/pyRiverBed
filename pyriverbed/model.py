"""The pyRiverBed model.

:class:`RiverBedModel` owns a :class:`~pyriverbed.config.Config` and runs the
workflow: build the planform, compute curvature, generate the bed, expand it to
a point cloud, write the files, and -- if migration is on -- step the planform
forward in time, detecting cutoffs as it goes.

All state lives on the instance, so several models can coexist in one process.
That was impossible in v1.x, where every parameter was a module global and the
Numba kernels baked them in at compile time.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

from . import bed as bed_module
from . import geometry, mesh, planform, plotting
from ._version import __version__
from .config import Config, default_config, load_config
from .logging_utils import (ProgressBar, format_duration, get_logger,
                            log_banner, log_key_values, log_table, stage)
from .migration import (Cutoff, carve_cutoff, find_chute_cutoff,
                        find_neck_cutoff, migrate)
from .planform import Centerline

log = get_logger(__name__)

__all__ = ["RiverBedModel", "RunResult", "run"]

#: Called as ``callback(step, n_steps)`` after each time step. Returning
#: ``False`` stops the run -- this is how the GUI implements its Stop button.
ProgressCallback = Callable[[int, int], bool | None]


@dataclass
class RunResult:
    """Everything a finished run produced.

    Attributes
    ----------
    config
        The configuration that was run.
    centerline
        The final centerline.
    bed
        The final bed topography on the ``(s, n)`` grid.
    cloud
        The final riverbed point cloud, ``(n, 3)``.
    curvature_original, curvature_filtered, curvature_lagged
        The final curvature signal at its three stages.
    sinuosity, migration_rate
        Time series, one value per migration time step.
    cutoffs
        Every cutoff that happened, in order.
    files
        Every file written.
    centerline_history
        Centerline coordinates at each graphic printout.
    steps_completed
        Time steps actually run, which is less than requested if stopped early.
    elapsed
        Wall-clock duration (s).
    """

    config: Config
    centerline: Centerline
    bed: bed_module.BedTopography
    cloud: np.ndarray
    curvature_original: np.ndarray
    curvature_filtered: np.ndarray
    curvature_lagged: np.ndarray
    sinuosity: np.ndarray = field(default_factory=lambda: np.zeros(0))
    migration_rate: np.ndarray = field(default_factory=lambda: np.zeros(0))
    cutoffs: list[Cutoff] = field(default_factory=list)
    files: list[Path] = field(default_factory=list)
    centerline_history: list[tuple[np.ndarray, np.ndarray]] = field(
        default_factory=list)
    steps_completed: int = 0
    elapsed: float = 0.0

    @property
    def neck_cutoffs(self) -> list[Cutoff]:
        """Only the neck cutoffs."""
        return [c for c in self.cutoffs if c.kind == "neck"]

    @property
    def chute_cutoffs(self) -> list[Cutoff]:
        """Only the chute cutoffs."""
        return [c for c in self.cutoffs if c.kind == "chute"]

    def summary(self) -> dict[str, object]:
        """A short dict of headline numbers, handy in a notebook."""
        return {
            "steps": self.steps_completed,
            "nodes": self.centerline.n_nodes,
            "length_m": round(self.centerline.length, 2),
            "sinuosity": round(self.centerline.sinuosity, 4),
            "bar_pool_relief_m": round(self.bed.relief, 4),
            "neck_cutoffs": len(self.neck_cutoffs),
            "chute_cutoffs": len(self.chute_cutoffs),
            "files": len(self.files),
            "elapsed_s": round(self.elapsed, 2),
        }


class RiverBedModel:
    """Generate synthetic riverbed topography for a meandering river.

    Parameters
    ----------
    config
        The configuration to run. ``None`` uses :func:`~pyriverbed.config.default_config`.

    Examples
    --------
    >>> from pyriverbed import RiverBedModel, default_config    # doctest: +SKIP
    >>> config = default_config()                              # doctest: +SKIP
    >>> config.channel.width = 0.8                             # doctest: +SKIP
    >>> result = RiverBedModel(config).run()                   # doctest: +SKIP
    >>> result.summary()['sinuosity']                          # doctest: +SKIP
    """

    def __init__(self, config: Config | None = None) -> None:
        self.config = (config or default_config()).validate()
        self.rng = np.random.default_rng(self.config.migration.seed)
        self.centerline: Centerline | None = None
        self.cutoffs: list[Cutoff] = []
        self.files: list[Path] = []
        self._frames: plotting.RunFrames | None = None
        self._history: list[tuple[np.ndarray, np.ndarray]] = []

    # -- classmethods -----------------------------------------------------

    @classmethod
    def from_file(cls, path: str | Path) -> "RiverBedModel":
        """Build a model from an input file, v2 or legacy v1.x."""
        return cls(load_config(path))

    # -- properties -------------------------------------------------------

    @property
    def output_dir(self) -> Path:
        """The output directory, created on first access."""
        path = self.config.output.path
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def stem(self) -> str:
        """Stem of every output file name."""
        return self.config.stem

    # -- reporting --------------------------------------------------------

    def log_configuration(self) -> None:
        """Log the parameter tables for this run."""
        config = self.config
        channel = config.channel

        if config.mode == "kinoshita":
            log.info("Mode: generate a Kinoshita curve from its equation")
            log.info("   %s", planform.kinoshita_equation_text(config.kinoshita))
        else:
            log.info("Mode: read a river centerline from %s",
                     config.centerline_file)

        rows = [
            ["Channel width", f"{channel.width:g}", "m"],
            ["Flow depth", f"{channel.depth:g}", "m"],
            ["Channel slope", f"{channel.slope:g}", "-"],
            ["Half-width / depth", f"{channel.beta:.3g}", "-"],
            ["Scour factor A", f"{channel.scour_factor:.4g}", "-"],
            ["Transverse slope corrector",
             f"{channel.transverse_slope_corrector:g}", "-"],
            ["Transverse resolution", f"{channel.interval:.4g}", "m"],
            ["Transverse # of nodes", f"{2 * channel.n_offsets + 1}", "-"],
            ["Curvature method", config.curvature.method, "-"],
            ["Smoothing passes (initial)", f"{config.curvature.n_passes}", "-"],
            ["Curvature phase lag",
             f"{config.lag.strength:g} widths" if config.lag.enabled else "off",
             "-"],
        ]
        if config.mode == "kinoshita":
            rows[7:7] = [["Streamwise resolution", f"{channel.ds:g}", "m"]]
        log_table(rows, headers=["Parameter", "Value", "Unit"],
                  align=("left", "right", "left"))

        if config.migration.enabled:
            m = config.migration
            log_table(
                [
                    ["Time steps", f"{m.n_steps}", "-"],
                    ["Time step dt", f"{m.dt:g}", "s"],
                    ["Simulated time", f"{m.n_steps * m.dt / 86400:.4g}", "d"],
                    ["Erosion coefficient E0", f"{m.e0:g}", "1/s"],
                    ["E0 * dt", f"{m.e0 * m.dt:g}", "widths/step"],
                    ["Inlet noise Ub0", f"{m.ub0:g}", "-"],
                    ["Inlet bias C0", f"{m.c0:g}", "-"],
                    ["Friction Cf0", f"{m.cf0:g}", "-"],
                    ["Froude Fr0", f"{m.fr0:g}", "-"],
                    ["Smoothing passes (per step)",
                     f"{config.curvature.n_migration_passes}", "-"],
                    ["Neck cutoffs",
                     "on" if config.neck_active else "off", "-"],
                    ["RNG seed",
                     "unseeded" if m.seed is None else str(m.seed), "-"],
                ],
                headers=["Migration parameter", "Value", "Unit"],
                align=("left", "right", "left"))

        if config.chute_active:
            c = config.chute_cutoff
            log_table(
                [
                    ["Frequency", f"{c.frequency * 100:g}%", "per time step"],
                    ["Mean recurrence",
                     f"{1 / c.frequency:.4g}" if c.frequency else "never",
                     "time steps"],
                    ["Starting step", f"{c.start_step}", "-"],
                    ["Entrance location", c.entrance, "-"],
                    ["Span", f"{c.span}", "# of entrances"],
                    ["Max angle with valley axis",
                     f"{c.max_valley_angle:g}", "deg"],
                    ["Min bypassed length",
                     f"{c.min_length_widths:g}", "channel widths"],
                    ["Min bypassed sinuosity", f"{c.min_sinuosity:g}", "-"],
                    ["Margin at reach ends", f"{c.end_margin}",
                     "# of entrances"],
                ],
                headers=["Chute cutoff parameter", "Value", "Unit"],
                align=("left", "right", "left"))

    def _log_planform(self, line: Centerline) -> None:
        """Log the geometry of the centerline that was built."""
        mean, median, mode = line.spacing_stats()
        log_key_values([
            ("nodes", line.n_nodes),
            ("channel length", f"{line.length:.4g} m"),
            ("valley length", f"{line.valley_length:.4g} m"),
            ("sinuosity", f"{line.sinuosity:.4f}"),
            ("node spacing", f"mean {mean:.4g} m, median {median:.4g} m, "
                             f"mode {mode:.4g} m"),
        ])

    # -- pipeline ---------------------------------------------------------

    def build_planform(self) -> Centerline:
        """Build the initial centerline and store it on the model."""
        with stage("Building channel planform"):
            line = planform.build_centerline(self.config)
        self._log_planform(line)
        self.centerline = line
        return line

    def compute_curvature_fields(
        self, line: Centerline
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the raw, filtered and phase-lagged curvature of *line*.

        The bed is built from the *lagged* curvature; the migration model uses
        both the filtered and the lagged one.
        """
        config = self.config
        cur_raw, theta = planform.curvature(line.s, line.x, line.y,
                                           config.curvature)
        line.curvature, line.theta = cur_raw, theta
        cur_filtered = planform.filter_curvature(cur_raw)
        if config.lag.enabled:
            spacing = line.mean_spacing
            window = (int(config.lag.strength * config.channel.width / spacing)
                      if spacing > 0 else 0)
            cur_lagged = planform.phase_lag(cur_filtered, window)
        else:
            cur_lagged = cur_filtered.copy()
        return cur_raw, cur_filtered, cur_lagged

    def build_bed(self, line: Centerline, cur_lagged: np.ndarray):
        """Build the bed topography and the point cloud for *line*."""
        bed = bed_module.compute_bed(cur_lagged, line.s, self.config.channel,
                                     self.config.flip)
        cloud, banklines = geometry.build_point_cloud(
            line.x, line.y, bed, self.config.channel,
            collect_banklines=self.config.output.save_bankline)
        return bed, cloud, banklines

    def write_outputs(self, line: Centerline, bed, cloud: np.ndarray,
                      banklines) -> list[Path]:
        """Write the point cloud, banklines and mesh files."""
        out = self.config.output
        directory = self.output_dir
        written: list[Path] = []
        if out.save_xyz:
            written.append(mesh.write_xyz(cloud, directory /
                                          f"{self.stem}_topo.xyz"))
        if out.save_bankline and banklines is not None:
            written.append(mesh.write_bankline(
                geometry.bankline_polygon(banklines),
                directory / f"{self.stem}_boundary.i2s"))
        if out.save_mesh:
            written += mesh.write_mesh_files(cloud, line.n_nodes,
                                            self.config.channel, directory,
                                            self.stem)
        self.files += written
        return written

    # -- figures ----------------------------------------------------------

    def _frame_dirs(self) -> tuple[Path, Path]:
        """Directories holding the two animation frame sequences."""
        base = self.output_dir
        return base / "frames_bed", base / "frames_planform"

    def _draw(self, line: Centerline, cloud: np.ndarray, bed,
              cur_raw: np.ndarray, cur_filtered: np.ndarray,
              cur_lagged: np.ndarray, step: int) -> None:
        """Draw and save the figures for one graphic printout."""
        out = self.config.output
        if not out.save_figures:
            return
        directory = self.output_dir
        migrating = self.config.migration.enabled

        headline: list[Path] = []
        if step == 0:
            headline = [directory / f"{self.stem}_pyriverbed.{ext}"
                        for ext in out.figure_formats]
        frame_paths = list(headline)
        bed_frame = planform_frame = None
        if migrating:
            bed_dir, plan_dir = self._frame_dirs()
            bed_frame = bed_dir / f"{self.stem}_{step:08d}_bed.png"
            planform_frame = plan_dir / f"{self.stem}_{step:08d}_planform.png"
            frame_paths.append(bed_frame)

        figure = plotting.plot_summary(
            line, cloud, bed, cur_raw, cur_filtered, cur_lagged, self.config,
            step, cutoffs=self.cutoffs, paths=frame_paths,
            dpi=out.figure_dpi if step == 0 and not migrating else out.frame_dpi)
        figure.clf()
        if headline:
            self.files += headline

        if migrating:
            figure = plotting.plot_meander_belt(
                line, self.config, step, cutoffs=self.cutoffs,
                paths=[planform_frame], dpi=out.frame_dpi,
                history_x=[h[0] for h in self._history],
                history_y=[h[1] for h in self._history])
            figure.clf()
            if self._frames is not None:
                self._frames.register(bed_frame, planform_frame)

    # -- the run ----------------------------------------------------------

    def run(self, progress: ProgressCallback | None = None) -> RunResult:
        """Execute the whole workflow.

        Parameters
        ----------
        progress
            Called as ``progress(step, n_steps)`` after each time step.
            Returning ``False`` stops the run early and still writes the
            diagnostics and animations for the steps completed.

        Returns
        -------
        RunResult
        """
        started = time.perf_counter()
        config = self.config
        n_steps = config.n_steps

        log.info("pyRiverBed %s starting in %s", __version__,
                 self.output_dir.resolve())
        self.log_configuration()

        line = self.build_planform()
        self._frames = plotting.RunFrames(directories=self._frame_dirs())
        sinuosity = np.zeros(n_steps)
        rate = np.zeros(n_steps)
        cur_raw = cur_filtered = cur_lagged = np.zeros(line.n_nodes)
        bed = cloud = None
        stopped_early = False

        bar = ProgressBar(max(n_steps, 1),
                          label="   migrating ") if n_steps else None
        for step in range(n_steps + 1):
            verbose = (step % config.migration.log_every == 0) or not n_steps
            cur_raw, cur_filtered, cur_lagged = self.compute_curvature_fields(line)
            bed, cloud, banklines = self.build_bed(line, cur_lagged)

            if step == 0:
                with stage("Writing riverbed, bankline and mesh files"):
                    self.write_outputs(line, bed, cloud, banklines)

            if step % config.migration.plot_every == 0 or step == n_steps:
                self._history.append((line.x.copy(), line.y.copy()))
                self._draw(line, cloud, bed, cur_raw, cur_filtered,
                           cur_lagged, step)

            if step == n_steps:
                break

            result = migrate(line, cur_filtered, cur_lagged, config.channel,
                             config.migration, self.rng)
            line = result.centerline
            rate[step] = result.mean_rate
            line = self._apply_cutoffs(line, cur_filtered, step)

            _, x, y = planform.smooth(line.x, line.y,
                                      config.curvature.n_migration_passes)
            s, x, y = planform.resample(x, y, config.channel.interval)
            cur, theta = planform.curvature(s, x, y, config.curvature)
            line = Centerline(x=x, y=y, s=s, curvature=cur, theta=theta)
            sinuosity[step] = line.sinuosity

            if verbose:
                log.log(15, "   step %d/%d  sinuosity %.4f  rate %.4g m/step "
                        " nodes %d  cutoffs %d", step + 1, n_steps,
                        line.sinuosity, rate[step], line.n_nodes,
                        len(self.cutoffs))
            if bar is not None:
                bar.update(step + 1,
                           f"sinuosity {line.sinuosity:.3f}, "
                           f"{len(self.cutoffs)} cutoffs")
            if progress is not None and progress(step + 1, n_steps) is False:
                log.warning("run stopped at step %d of %d on request",
                            step + 1, n_steps)
                stopped_early = True
                sinuosity = sinuosity[:step + 1]
                rate = rate[:step + 1]
                break
        if bar is not None:
            bar.close()

        self.centerline = line
        steps_done = len(sinuosity) if stopped_early else n_steps
        self._finalise(sinuosity, rate)
        elapsed = time.perf_counter() - started

        result = RunResult(
            config=config, centerline=line, bed=bed, cloud=cloud,
            curvature_original=cur_raw, curvature_filtered=cur_filtered,
            curvature_lagged=cur_lagged, sinuosity=sinuosity,
            migration_rate=rate, cutoffs=list(self.cutoffs),
            files=list(self.files), centerline_history=list(self._history),
            steps_completed=steps_done, elapsed=elapsed)
        self._log_summary(result)
        return result

    def _apply_cutoffs(self, line: Centerline, cur_filtered: np.ndarray,
                       step: int) -> Centerline:
        """Detect and carve cutoffs for one time step.

        A neck cutoff is looked for first. A chute cutoff is only looked for
        when no neck cutoff fired, because carving one invalidates the node
        indices the other's search returned.
        """
        config = self.config
        i, j = find_neck_cutoff(line, config.channel, config.neck_cutoff)
        if i >= 0:
            line, cut = carve_cutoff(line, i, j, "neck", step, config.channel)
            self.cutoffs.append(cut)
            log.info("   neck cutoff at step %d: nodes %d-%d bypassed "
                     "(%d-node oxbow lake)", step, i, j, cut.bypassed_nodes)
            return line

        if config.chute_active:
            i, j = find_chute_cutoff(line, cur_filtered, step, config.channel,
                                     config.chute_cutoff, self.rng)
            if i >= 0:
                line, cut = carve_cutoff(line, i, j, "chute", step, config.channel)
                self.cutoffs.append(cut)
                log.info("   chute cutoff at step %d: nodes %d-%d bypassed "
                         "(%d-node oxbow lake)", step, i, j,
                         cut.bypassed_nodes)
        return line

    def _finalise(self, sinuosity: np.ndarray, rate: np.ndarray) -> None:
        """Write the time series, the diagnostics figure and the animations."""
        config = self.config
        if not config.migration.enabled or sinuosity.size == 0:
            return
        directory = self.output_dir
        for name, data in (("sinuosity", sinuosity),
                           ("mean_migration_rate", rate)):
            path = directory / f"{self.stem}_{name}.txt"
            np.savetxt(path, data, fmt="%.8e")
            self.files.append(path)
        log.info("   wrote %s_sinuosity.txt, %s_mean_migration_rate.txt",
                 self.stem, self.stem)

        if config.output.save_figures:
            paths = [directory / f"{self.stem}_diagnostics.png"]
            figure = plotting.plot_diagnostics(sinuosity, rate, self.cutoffs,
                                               config, paths=paths,
                                               dpi=config.output.figure_dpi)
            figure.clf()
            self.files += paths
            log.info("   wrote %s", paths[0].name)

        if config.output.save_gif and self._frames is not None:
            with stage("Assembling animations"):
                self.files += plotting.make_gifs(
                    self._frames, directory, self.stem,
                    fps=config.output.gif_fps)

    def _log_summary(self, result: RunResult) -> None:
        """Log the closing summary of a run."""
        log.info("")
        rows = [[key.replace("_", " "), str(value)]
                for key, value in result.summary().items()]
        log_table(rows, headers=["Result", "Value"], align=("left", "right"))
        if result.cutoffs:
            log.info("   %d neck and %d chute cutoff(s)",
                     len(result.neck_cutoffs), len(result.chute_cutoffs))
        log.info("Finished in %s", format_duration(result.elapsed))


def run(config: Config | str | Path | None = None,
        progress: ProgressCallback | None = None) -> RunResult:
    """Run pyRiverBed.

    The one-call entry point used by the CLI, the GUI and the notebook helpers.

    Parameters
    ----------
    config
        A :class:`~pyriverbed.config.Config`, a path to an input file, or
        ``None`` for the defaults.
    progress
        Progress callback, see :meth:`RiverBedModel.run`.

    Returns
    -------
    RunResult
    """
    if config is None or isinstance(config, Config):
        model = RiverBedModel(config)
    else:
        model = RiverBedModel.from_file(config)
    return model.run(progress=progress)
