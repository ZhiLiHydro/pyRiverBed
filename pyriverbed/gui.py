"""Graphical interface for pyRiverBed.

Still tkinter -- it ships with Python, and adding Qt would undo the work of
getting the dependency list down to three packages -- but rebuilt around
``ttk``, a tabbed notebook and a grid layout instead of the v1.x wall of
hand-packed frames.

What is new beyond the looks:

* the model runs on a **worker thread**, so the window stays responsive and a
  long migration can be stopped;
* log records are piped into a **live log pane** through a queue, so the run
  narrates itself in the window instead of in a terminal behind it;
* a **live preview** of the planform, drawn with the Matplotlib that is already
  a dependency;
* every field is **validated as you type**, and input files can be loaded and
  saved.
"""

from __future__ import annotations

import logging
import queue
import threading
import traceback
import webbrowser
from dataclasses import fields, is_dataclass
from pathlib import Path
from typing import Any, Callable

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from ._version import __version__
from .config import (ChuteCutoffConfig, Config, ConfigError, CurvatureConfig,
                     default_config, load_config, write_config)
from .logging_utils import LOGGER_NAME, TRACE, get_logger

log = get_logger(__name__)

REPO_URL = "https://github.com/ZhiLiHydro/pyRiverBed"

#: Field spec: (section, key, label, widget kind, choices/limits)
_BOOL = "bool"
_NUM = "num"
_INT = "int"
_TEXT = "text"
_CHOICE = "choice"
_FILE = "file"
_DIR = "dir"
_SCALE = "scale"


class _LogHandler(logging.Handler):
    """Push formatted log records onto a queue for the GUI to drain."""

    def __init__(self, sink: queue.Queue) -> None:
        super().__init__()
        self.sink = sink

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self.sink.put_nowait((record.levelno, self.format(record)))
        except queue.Full:  # pragma: no cover - the queue is unbounded
            pass


class PyRiverBedGUI:
    """The pyRiverBed main window.

    Parameters
    ----------
    root
        The Tk root window.
    config
        Configuration to pre-load. ``None`` uses the defaults, which reproduce
        the Abad & Garcia (2009) laboratory flume.
    """

    def __init__(self, root: tk.Tk, config: Config | None = None) -> None:
        self.root = root
        self.config = config or default_config()
        self.vars: dict[tuple[str, str], tk.Variable] = {}
        self.errors: dict[tuple[str, str], tk.Label] = {}
        self.scale_labels: dict[tuple[str, str], tk.Label] = {}
        self.log_queue: queue.Queue = queue.Queue()
        self.worker: threading.Thread | None = None
        self.stop_requested = threading.Event()
        self.input_path: Path | None = None
        self._preview_canvas = None
        self._preview_figure = None

        root.title(f"pyRiverBed {__version__}")
        root.minsize(940, 700)
        self._init_style()
        self._build()
        self._config_to_widgets(self.config)
        self._attach_log_handler()
        self._drain_log()
        root.protocol("WM_DELETE_WINDOW", self._on_close)
        log.info("pyRiverBed %s ready. Set the parameters, then press "
                 "'Run pyRiverBed'.", __version__)
        log.info("Press 'Preview planform' to check the centerline first.")

    # -- chrome -----------------------------------------------------------

    def _init_style(self) -> None:
        """Pick the best available ttk theme and set up named styles.

        Font sizes are derived from the platform default rather than hard-coded,
        so the window looks right whatever the system font happens to be.
        """
        import tkinter.font as tkfont

        style = ttk.Style(self.root)
        for candidate in ("clam", "alt", "vista", "aqua", "default"):
            if candidate in style.theme_names():
                style.theme_use(candidate)
                break
        base = tkfont.nametofont("TkDefaultFont")
        family = base.actual("family")
        size = abs(base.actual("size")) or 10
        style.configure("Header.TLabel", font=(family, size + 7, "bold"))
        style.configure("Sub.TLabel", foreground="#555555",
                        font=(family, size))
        style.configure("Error.TLabel", foreground="#b00020",
                        font=(family, max(size - 1, 7)))
        style.configure("Run.TButton", font=(family, size + 1, "bold"))
        style.configure("Group.TLabelframe.Label", font=(family, size, "bold"))

    def _build(self) -> None:
        """Lay out the whole window."""
        root = self.root
        root.columnconfigure(0, weight=1)
        root.rowconfigure(1, weight=1)

        header = ttk.Frame(root, padding=(12, 10, 12, 6))
        header.grid(row=0, column=0, sticky="ew")
        header.columnconfigure(1, weight=1)
        ttk.Label(header, text="pyRiverBed", style="Header.TLabel").grid(
            row=0, column=0, sticky="w")
        ttk.Label(header,
                  text="Synthetic riverbed topography for meandering rivers",
                  style="Sub.TLabel").grid(row=1, column=0, sticky="w")
        buttons = ttk.Frame(header)
        buttons.grid(row=0, column=2, rowspan=2, sticky="e")
        ttk.Button(buttons, text="Open...", command=self.on_open).pack(
            side="left", padx=2)
        ttk.Button(buttons, text="Save as...", command=self.on_save).pack(
            side="left", padx=2)
        ttk.Button(buttons, text="Defaults", command=self.on_defaults).pack(
            side="left", padx=2)
        ttk.Button(buttons, text="GitHub", command=self.on_github).pack(
            side="left", padx=2)

        panes = ttk.Panedwindow(root, orient="horizontal")
        panes.grid(row=1, column=0, sticky="nsew", padx=10, pady=4)

        left = ttk.Frame(panes)
        self.notebook = ttk.Notebook(left)
        self.notebook.pack(fill="both", expand=True)
        self._build_tabs()
        panes.add(left, weight=3)

        right = ttk.Frame(panes)
        right.rowconfigure(0, weight=3)
        right.rowconfigure(1, weight=2)
        right.columnconfigure(0, weight=1)
        self._build_preview(right)
        self._build_log(right)
        panes.add(right, weight=4)

        footer = ttk.Frame(root, padding=(12, 4, 12, 10))
        footer.grid(row=2, column=0, sticky="ew")
        footer.columnconfigure(1, weight=1)
        self.run_button = ttk.Button(footer, text="Run pyRiverBed",
                                     style="Run.TButton", command=self.on_run)
        self.run_button.grid(row=0, column=0, sticky="w")
        self.stop_button = ttk.Button(footer, text="Stop", command=self.on_stop,
                                      state="disabled")
        self.stop_button.grid(row=0, column=1, sticky="w", padx=6)
        self.progress = ttk.Progressbar(footer, mode="determinate")
        self.progress.grid(row=0, column=2, sticky="ew", padx=10)
        footer.columnconfigure(2, weight=1)
        self.status = ttk.Label(footer, text="ready", style="Sub.TLabel")
        self.status.grid(row=0, column=3, sticky="e")

    def _build_preview(self, parent) -> None:
        """Create the planform preview panel."""
        box = ttk.Labelframe(parent, text="Preview", padding=6,
                             style="Group.TLabelframe")
        box.grid(row=0, column=0, sticky="nsew", pady=(0, 6))
        box.rowconfigure(0, weight=1)
        box.columnconfigure(0, weight=1)
        try:
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            from matplotlib.figure import Figure
        except ImportError:  # pragma: no cover
            ttk.Label(box, text="matplotlib is unavailable, so no preview",
                      style="Sub.TLabel").grid(row=0, column=0)
            return
        figure = Figure(figsize=(5, 3.2), layout="constrained")
        self._preview_axes = figure.add_subplot(111)
        self._preview_axes.set_title("press Preview to draw the planform",
                                     fontsize=9)
        self._preview_axes.tick_params(labelsize=8)
        canvas = FigureCanvasTkAgg(figure, master=box)
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self._preview_figure = figure
        self._preview_canvas = canvas
        ttk.Button(box, text="Preview planform",
                   command=self.on_preview).grid(row=1, column=0, sticky="w",
                                                 pady=(6, 0))

    def _build_log(self, parent) -> None:
        """Create the live log pane."""
        box = ttk.Labelframe(parent, text="Log", padding=6,
                             style="Group.TLabelframe")
        box.grid(row=1, column=0, sticky="nsew")
        box.rowconfigure(0, weight=1)
        box.columnconfigure(0, weight=1)
        self.log_text = tk.Text(box, height=12, wrap="none", state="disabled",
                                font=("TkFixedFont", 9), background="#1e1e1e",
                                foreground="#d4d4d4", insertbackground="#d4d4d4",
                                borderwidth=0)
        self.log_text.grid(row=0, column=0, sticky="nsew")
        yscroll = ttk.Scrollbar(box, orient="vertical",
                                command=self.log_text.yview)
        yscroll.grid(row=0, column=1, sticky="ns")
        xscroll = ttk.Scrollbar(box, orient="horizontal",
                                command=self.log_text.xview)
        xscroll.grid(row=1, column=0, sticky="ew")
        self.log_text.configure(yscrollcommand=yscroll.set,
                                xscrollcommand=xscroll.set)
        for tag, colour in (("debug", "#808080"), ("info", "#d4d4d4"),
                            ("warning", "#dcdcaa"), ("error", "#f48771")):
            self.log_text.tag_configure(tag, foreground=colour)
        ttk.Button(box, text="Clear", command=self.on_clear_log).grid(
            row=2, column=0, sticky="w", pady=(6, 0))

    # -- field construction ----------------------------------------------

    def _build_tabs(self) -> None:
        """Populate the notebook with one tab per group of parameters."""
        self._tab("Planform", [
            ("Mode", [
                ("", "mode", "Mode", _CHOICE, Config.MODES),
                ("", "centerline_file", "Centerline file", _FILE, None),
                ("", "name", "Output name (optional)", _TEXT, None),
                ("curvature", "smoothing_level", "Smoothing level", _SCALE,
                 (0, 100)),
            ]),
            ("Kinoshita curve (mode = kinoshita)", [
                ("kinoshita", "n_bends", "Number of bends", _INT, None),
                ("kinoshita", "arc_wavelength", "Arc wavelength (m)", _NUM, None),
                ("kinoshita", "max_angular_amplitude",
                 "Max angular amplitude (deg)", _NUM, None),
                ("kinoshita", "skewness", "Skewness", _NUM, None),
                ("kinoshita", "flatness", "Flatness", _NUM, None),
            ]),
        ])
        self._tab("Channel & bed", [
            ("Channel geometry", [
                ("channel", "width", "Width (m)", _NUM, None),
                ("channel", "depth", "Depth (m)", _NUM, None),
                ("channel", "slope", "Slope", _NUM, None),
                ("channel", "transverse_slope_corrector",
                 "Transverse slope corrector", _NUM, None),
            ]),
            ("Resolution", [
                ("channel", "ds", "Streamwise spacing (m)", _NUM, None),
                ("channel", "n_offsets", "Offsets per side", _INT, None),
            ]),
            ("Curvature & phase lag", [
                ("curvature", "method", "Curvature method", _CHOICE,
                 CurvatureConfig.METHODS),
                ("curvature", "despike", "De-spike curvature", _BOOL, None),
                ("lag", "enabled", "Curvature phase lag", _BOOL, None),
                ("lag", "strength", "Lag strength (channel widths)", _NUM, None),
            ]),
            ("Flip", [
                ("flip", "streamwise", "Flip streamwise", _BOOL, None),
                ("flip", "transverse", "Flip transverse", _BOOL, None),
            ]),
        ])
        self._tab("Migration", [
            ("Migration", [
                ("migration", "enabled", "Meander migration", _BOOL, None),
                ("migration", "n_steps", "Number of time steps", _INT, None),
                ("migration", "dt", "Time step dt (s)", _NUM, None),
                ("migration", "e0", "Erosion coefficient E0 (1/s)", _NUM, None),
                ("migration", "seed", "RNG seed (blank = unseeded)", _TEXT,
                 None),
            ]),
            ("Hydraulics", [
                ("migration", "ub0", "Inlet noise Ub0", _NUM, None),
                ("migration", "c0", "Inlet bias C0", _NUM, None),
                ("migration", "cf0", "Friction Cf0", _NUM, None),
                ("migration", "fr0", "Froude Fr0", _NUM, None),
            ]),
            ("Stability", [
                ("curvature", "migration_smoothing_level",
                 "Smoothing per step (2-8)", _INT, None),
                ("migration", "end_taper_widths",
                 "Fixed-end taper (widths)", _NUM, None),
            ]),
            ("Reporting", [
                ("migration", "log_every", "Log every N steps", _INT, None),
                ("migration", "plot_every", "Plot every N steps", _INT, None),
            ]),
        ])
        self._tab("Cutoffs", [
            ("Neck cutoff (geometric, deterministic)", [
                ("neck_cutoff", "enabled", "Detect neck cutoffs", _BOOL, None),
                ("neck_cutoff", "min_separation_widths",
                 "Min along-channel separation (widths)", _NUM, None),
                ("neck_cutoff", "max_distance_widths",
                 "Max node distance (widths)", _NUM, None),
                ("neck_cutoff", "end_margin_widths",
                 "Margin at reach ends (widths)", _NUM, None),
            ]),
            ("Chute cutoff (geometric criteria + random trigger)", [
                ("chute_cutoff", "enabled", "Model chute cutoffs", _BOOL, None),
                ("chute_cutoff", "frequency",
                 "Frequency (probability per step)", _NUM, None),
                ("chute_cutoff", "start_step", "Starting time step", _INT, None),
                ("chute_cutoff", "entrance", "Entrance location", _CHOICE,
                 ChuteCutoffConfig.ENTRANCES),
                ("chute_cutoff", "span", "Span (# of entrances)", _INT, None),
                ("chute_cutoff", "max_valley_angle",
                 "Max angle with valley axis (deg)", _NUM, None),
                ("chute_cutoff", "min_length_widths",
                 "Min bypassed length (widths)", _NUM, None),
                ("chute_cutoff", "min_sinuosity",
                 "Min bypassed sinuosity", _NUM, None),
                ("chute_cutoff", "end_margin", "Margin at reach ends", _INT,
                 None),
            ]),
        ])
        self._tab("Output", [
            ("Where", [
                ("output", "directory", "Output directory", _DIR, None),
                ("output", "log_file", "Log file name", _TEXT, None),
            ]),
            ("What", [
                ("output", "save_xyz", "Point cloud (.xyz)", _BOOL, None),
                ("output", "save_bankline", "Banklines (.i2s)", _BOOL, None),
                ("output", "save_mesh", "FEM mesh & BC files", _BOOL, None),
                ("output", "save_figures", "Figures", _BOOL, None),
                ("output", "save_gif", "Animated GIFs", _BOOL, None),
            ]),
            ("Figure settings", [
                ("output", "figure_formats", "Headline figure formats", _TEXT,
                 None),
                ("output", "figure_dpi", "Figure DPI", _INT, None),
                ("output", "frame_dpi", "Animation frame DPI", _INT, None),
                ("output", "gif_fps", "GIF frames per second", _INT, None),
            ]),
        ])

    def _tab(self, title: str, groups) -> None:
        """Add one notebook tab holding *groups* of fields."""
        outer = ttk.Frame(self.notebook)
        self.notebook.add(outer, text=title)
        canvas = tk.Canvas(outer, borderwidth=0, highlightthickness=0)
        scroll = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        inner = ttk.Frame(canvas, padding=8)
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")
        window = canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>",
                   lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>",
                    lambda e: canvas.itemconfigure(window, width=e.width))

        for group_title, specs in groups:
            box = ttk.Labelframe(inner, text=group_title, padding=8,
                                 style="Group.TLabelframe")
            box.pack(fill="x", pady=(0, 8))
            box.columnconfigure(1, weight=1)
            for row, spec in enumerate(specs):
                self._field(box, row, *spec)

    def _field(self, parent, row: int, section: str, key: str, label: str,
               kind: str, extra: Any) -> None:
        """Create one labelled input widget."""
        ident = (section, key)
        ttk.Label(parent, text=label).grid(row=row * 2, column=0, sticky="w",
                                          padx=(0, 8), pady=2)
        if kind == _BOOL:
            var: tk.Variable = tk.BooleanVar()
            widget = ttk.Checkbutton(parent, variable=var, text="")
        elif kind == _CHOICE:
            var = tk.StringVar()
            widget = ttk.Combobox(parent, textvariable=var, values=list(extra),
                                  state="readonly", width=18)
        elif kind == _SCALE:
            var = tk.IntVar()
            holder = ttk.Frame(parent)
            value_label = ttk.Label(holder, width=4, style="Sub.TLabel")
            scale = ttk.Scale(holder, from_=extra[0], to=extra[1],
                              orient="horizontal", variable=var,
                              command=lambda v, lbl=value_label:
                              lbl.configure(text=str(int(float(v)))))
            scale.pack(side="left", fill="x", expand=True)
            value_label.pack(side="left", padx=(6, 0))
            self.scale_labels[ident] = value_label
            widget = holder
        elif kind in (_FILE, _DIR):
            var = tk.StringVar()
            holder = ttk.Frame(parent)
            entry = ttk.Entry(holder, textvariable=var)
            entry.pack(side="left", fill="x", expand=True)
            command = (self._pick_file if kind == _FILE else self._pick_dir)
            ttk.Button(holder, text="...", width=3,
                       command=lambda v=var, c=command: c(v)).pack(
                           side="left", padx=(4, 0))
            widget = holder
        else:
            var = tk.StringVar()
            widget = ttk.Entry(parent, textvariable=var)
        widget.grid(row=row * 2, column=1, sticky="ew", pady=1)
        self.vars[ident] = var

        # The error label is laid out but hidden, so a valid form has no gaps
        # and a message appears in place the moment a field goes bad.
        error = ttk.Label(parent, text="", style="Error.TLabel")
        error.grid(row=row * 2 + 1, column=1, sticky="w")
        error.grid_remove()
        self.errors[ident] = error
        if kind in (_NUM, _INT):
            var.trace_add("write",
                          lambda *_, i=ident, k=kind: self._check(i, k))

    def _pick_file(self, var: tk.Variable) -> None:
        """Ask for a centerline file."""
        path = filedialog.askopenfilename(
            title="Select a centerline file",
            filetypes=[("Text files", "*.txt *.dat *.csv"),
                       ("All files", "*.*")])
        if path:
            var.set(path)

    def _pick_dir(self, var: tk.Variable) -> None:
        """Ask for an output directory."""
        path = filedialog.askdirectory(title="Select an output directory")
        if path:
            var.set(path)

    def _check(self, ident: tuple[str, str], kind: str) -> bool:
        """Validate one numeric field and show or clear its error label."""
        text = str(self.vars[ident].get()).strip()
        label = self.errors[ident]
        if text == "":
            label.grid_remove()
            return True
        try:
            value = float(text)
            if kind == _INT and value != int(value):
                raise ValueError
        except ValueError:
            label.configure(text="not a number" if kind == _NUM
                            else "not a whole number")
            label.grid()
            return False
        label.grid_remove()
        return True

    # -- config <-> widgets ----------------------------------------------

    def _config_to_widgets(self, config: Config) -> None:
        """Copy *config* into the widgets."""
        data = config.to_dict()
        for (section, key), var in self.vars.items():
            value = data[key] if section == "" else data[section][key]
            if isinstance(var, tk.BooleanVar):
                var.set(bool(value))
            elif value is None:
                var.set("")
            elif isinstance(value, (tuple, list)):
                var.set(", ".join(str(v) for v in value))
            elif isinstance(var, tk.IntVar):
                var.set(int(value))
            else:
                var.set(f"{value:g}" if isinstance(value, float)
                        else str(value))
        for ident, label in self.scale_labels.items():
            label.configure(text=str(int(float(self.vars[ident].get()))))

    def _widgets_to_config(self) -> Config:
        """Build a validated :class:`Config` from the widgets.

        Raises
        ------
        ConfigError
            If any field cannot be read or is out of range.
        """
        data: dict[str, Any] = {}
        for (section, key), var in self.vars.items():
            value = var.get()
            if isinstance(value, str):
                value = value.strip()
            if section == "":
                data[key] = value
            else:
                data.setdefault(section, {})[key] = value
        # A blank seed means "unseeded", and a blank name means "derive it".
        if data.get("migration", {}).get("seed", "") == "":
            data["migration"]["seed"] = None
        return Config.from_dict(data).validate()

    # -- logging ----------------------------------------------------------

    def _attach_log_handler(self) -> None:
        """Route the pyRiverBed logger into the log pane."""
        handler = _LogHandler(self.log_queue)
        handler.setFormatter(logging.Formatter("%(message)s"))
        handler.setLevel(TRACE)
        logger = logging.getLogger(LOGGER_NAME)
        logger.setLevel(TRACE)
        logger.addHandler(handler)
        self._handler = handler

    def _drain_log(self) -> None:
        """Move queued log records into the text pane; reschedule itself."""
        pending = []
        try:
            while True:
                pending.append(self.log_queue.get_nowait())
        except queue.Empty:
            pass
        if pending:
            self.log_text.configure(state="normal")
            for level, message in pending:
                if level >= logging.ERROR:
                    tag = "error"
                elif level >= logging.WARNING:
                    tag = "warning"
                elif level >= TRACE:
                    tag = "info"
                else:
                    tag = "debug"
                self.log_text.insert("end", message + "\n", tag)
            self.log_text.see("end")
            self.log_text.configure(state="disabled")
        self.root.after(120, self._drain_log)

    def on_clear_log(self) -> None:
        """Empty the log pane."""
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.configure(state="disabled")

    # -- commands ---------------------------------------------------------

    def on_open(self) -> None:
        """Load an input file into the widgets."""
        path = filedialog.askopenfilename(
            title="Open a pyRiverBed input file",
            filetypes=[("pyRiverBed input", "*.ini *.cfg *.conf"),
                       ("v1.x steering file", "steering.txt *.txt"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            config = load_config(path)
        except (ConfigError, OSError, ValueError) as exc:
            messagebox.showerror("Cannot open", str(exc))
            return
        self.config = config
        self.input_path = Path(path)
        self._config_to_widgets(config)
        self.status.configure(text=f"loaded {Path(path).name}")
        log.info("loaded configuration from %s", path)

    def on_save(self) -> None:
        """Write the widgets out as an input file."""
        try:
            config = self._widgets_to_config()
        except ConfigError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return
        path = filedialog.asksaveasfilename(
            title="Save pyRiverBed input file", defaultextension=".ini",
            initialfile="pyriverbed.ini",
            filetypes=[("pyRiverBed input", "*.ini"), ("All files", "*.*")])
        if not path:
            return
        write_config(config, path)
        self.input_path = Path(path)
        self.status.configure(text=f"saved {Path(path).name}")

    def on_defaults(self) -> None:
        """Reset every field to the defaults."""
        if messagebox.askyesno("Reset", "Reset every parameter to its default?"):
            self.config = default_config()
            self._config_to_widgets(self.config)
            self.status.configure(text="defaults restored")

    def on_github(self) -> None:
        """Open the project page in a browser."""
        webbrowser.open(REPO_URL)

    def on_preview(self) -> None:
        """Draw the planform for the current settings, without running."""
        if self._preview_canvas is None:
            return
        try:
            config = self._widgets_to_config()
            from .planform import build_centerline
            line = build_centerline(config)
        except (ConfigError, OSError, ValueError) as exc:
            messagebox.showerror("Cannot preview", str(exc))
            return
        axes = self._preview_axes
        axes.clear()
        axes.plot(line.x, line.y, color="#c0392b", linewidth=1)
        axes.set_aspect("equal", adjustable="datalim")
        axes.set_title(f"{line.n_nodes} nodes, length {line.length:.3g} m, "
                       f"sinuosity {line.sinuosity:.3f}", fontsize=9)
        axes.tick_params(labelsize=8)
        self._preview_canvas.draw_idle()
        self.status.configure(text="preview drawn")

    def on_run(self) -> None:
        """Start a run on a worker thread."""
        if self.worker is not None and self.worker.is_alive():
            messagebox.showinfo("Already running",
                                "A run is already in progress.")
            return
        try:
            config = self._widgets_to_config()
        except ConfigError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return
        self.config = config
        self.stop_requested.clear()
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        total = max(config.n_steps, 1)
        self.progress.configure(maximum=total, value=0)
        self.status.configure(text="running...")
        self.worker = threading.Thread(target=self._worker, args=(config,),
                                       daemon=True)
        self.worker.start()

    def on_stop(self) -> None:
        """Ask the running model to stop at the next time step."""
        self.stop_requested.set()
        self.status.configure(text="stopping...")

    def _worker(self, config: Config) -> None:
        """Run the model off the UI thread."""
        from .model import RiverBedModel

        def progress(step: int, total: int) -> bool:
            self.root.after(0, self._set_progress, step, total)
            return not self.stop_requested.is_set()

        try:
            result = RiverBedModel(config).run(progress=progress)
        except Exception as exc:                          # noqa: BLE001
            log.error("run failed: %s", exc)
            log.debug("%s", traceback.format_exc())
            self.root.after(0, self._finish, None, exc)
            return
        self.root.after(0, self._finish, result, None)

    def _set_progress(self, step: int, total: int) -> None:
        """Update the progress bar from the UI thread."""
        self.progress.configure(maximum=max(total, 1), value=step)
        self.status.configure(text=f"step {step}/{total}")

    def _finish(self, result, error) -> None:
        """Re-enable the controls once a run ends."""
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        if error is not None:
            self.status.configure(text="failed")
            messagebox.showerror("Run failed", str(error))
            return
        self.progress.configure(value=self.progress["maximum"])
        summary = result.summary()
        self.status.configure(text=f"done in {summary['elapsed_s']} s")
        self._show_result(result)
        messagebox.showinfo(
            "Finished",
            "\n".join(f"{k.replace('_', ' ')}: {v}"
                      for k, v in summary.items())
            + f"\n\nOutput directory:\n{result.config.output.path.resolve()}")

    def _show_result(self, result) -> None:
        """Draw the final planform in the preview panel."""
        if self._preview_canvas is None:
            return
        axes = self._preview_axes
        axes.clear()
        for cut in result.cutoffs:
            if cut.oxbow_x.size:
                axes.fill(cut.oxbow_x, cut.oxbow_y, color="tan", zorder=0)
        line = result.centerline
        axes.plot(line.x, line.y, color="#2471a3", linewidth=1.2, zorder=1)
        axes.set_aspect("equal", adjustable="datalim")
        axes.set_title(f"final planform: sinuosity {line.sinuosity:.3f}, "
                       f"{len(result.cutoffs)} cutoff(s)", fontsize=9)
        axes.tick_params(labelsize=8)
        self._preview_canvas.draw_idle()

    def _on_close(self) -> None:
        """Confirm before closing while a run is in progress."""
        if self.worker is not None and self.worker.is_alive():
            if not messagebox.askyesno("Quit",
                                       "A run is in progress. Stop it and "
                                       "quit?"):
                return
            self.stop_requested.set()
        logging.getLogger(LOGGER_NAME).removeHandler(self._handler)
        self.root.destroy()


def main(config: Config | None = None) -> int:
    """Open the pyRiverBed GUI.

    Parameters
    ----------
    config
        Configuration to pre-load.

    Returns
    -------
    int
        Process exit status.
    """
    root = tk.Tk()
    PyRiverBedGUI(root, config)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
