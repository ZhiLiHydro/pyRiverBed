"""Screen output and logging for pyRiverBed.

Everything pyRiverBed says goes through :mod:`logging`, so a run can be
followed on screen, captured to a file, piped into a GUI pane or silenced
entirely without any of the computational code knowing about it.

The module also holds the small presentation helpers -- banner, tables,
progress bar, stage timer -- that make a run readable. They deliberately use
nothing but the standard library.
"""

from __future__ import annotations

import logging
import os
import shutil
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Sequence

from ._version import __version__

LOGGER_NAME = "pyriverbed"

#: Extra log level for the per-time-step chatter of a long migration run.
#: Sits between DEBUG and INFO so that ``-v`` can turn it on without also
#: turning on the very verbose DEBUG output.
TRACE = 15
logging.addLevelName(TRACE, "TRACE")


def get_logger(name: str | None = None) -> logging.Logger:
    """Return the pyRiverBed logger, or one of its children.

    Parameters
    ----------
    name
        Child name, typically ``__name__`` of the calling module. ``None``
        returns the root pyRiverBed logger.
    """
    if name is None or name == LOGGER_NAME:
        return logging.getLogger(LOGGER_NAME)
    return logging.getLogger(LOGGER_NAME).getChild(name.rsplit(".", 1)[-1])


# --------------------------------------------------------------------------
# colours
# --------------------------------------------------------------------------

class _Ansi:
    """ANSI escapes, blanked out when the stream cannot show them."""

    RESET = "\033[0m"
    DIM = "\033[2m"
    BOLD = "\033[1m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    CYAN = "\033[36m"

    def __init__(self, enabled: bool) -> None:
        if not enabled:
            for key in dir(self):
                if key.isupper():
                    setattr(self, key, "")


def _colour_supported(stream) -> bool:
    """Guess whether *stream* can render ANSI colour."""
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("PYRIVERBED_FORCE_COLOR"):
        return True
    if not hasattr(stream, "isatty") or not stream.isatty():
        return False
    if sys.platform == "win32":
        # Modern Windows terminals do, the legacy console does not. Enabling
        # virtual terminal processing is cheap, so try it and believe it.
        try:
            import ctypes

            kernel32 = ctypes.windll.kernel32
            return bool(kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7))
        except Exception:
            return False
    return os.environ.get("TERM", "") not in ("", "dumb")


class ConsoleFormatter(logging.Formatter):
    """Compact, coloured formatter for interactive use.

    ``INFO`` and ``TRACE`` records -- the normal narration of a run -- are
    printed as-is, so the output reads as prose rather than as a log. Anything
    unusual gets a visible, coloured tag.
    """

    def __init__(self, colour: bool = True) -> None:
        super().__init__()
        self.c = _Ansi(colour)

    def format(self, record: logging.LogRecord) -> str:
        message = record.getMessage()
        c = self.c
        if record.levelno >= logging.CRITICAL:
            head = f"{c.BOLD}{c.RED}[fatal]{c.RESET} "
        elif record.levelno >= logging.ERROR:
            head = f"{c.RED}[error]{c.RESET} "
        elif record.levelno >= logging.WARNING:
            head = f"{c.YELLOW}[warn]{c.RESET} "
        elif record.levelno >= logging.INFO:
            head = ""
        elif record.levelno >= TRACE:
            head = ""
        else:
            head = f"{c.DIM}[debug]{c.RESET} "
        text = head + message
        if record.exc_info:
            text += "\n" + self.formatException(record.exc_info)
        return text


class FileFormatter(logging.Formatter):
    """Timestamped, level-tagged formatter for the log file."""

    def __init__(self) -> None:
        super().__init__(
            fmt="%(asctime)s  %(levelname)-8s %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )


def setup_logging(
    level: int = logging.INFO,
    log_file: str | os.PathLike[str] | None = None,
    colour: bool | None = None,
    stream=None,
) -> logging.Logger:
    """Configure the pyRiverBed logger.

    Safe to call more than once: previous pyRiverBed handlers are removed
    first, so a notebook cell can be re-run without doubling every message.

    Parameters
    ----------
    level
        Threshold for the console. The log file, if any, always records
        everything from ``TRACE`` up.
    log_file
        Path of a log file to write in addition to the console. ``None``
        disables file logging.
    colour
        Force ANSI colour on or off. ``None`` auto-detects.
    stream
        Console stream. Defaults to :data:`sys.stderr`.

    Returns
    -------
    logging.Logger
        The configured pyRiverBed logger.
    """
    logger = logging.getLogger(LOGGER_NAME)
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    logger.setLevel(TRACE if level > TRACE else level)
    logger.propagate = False

    stream = sys.stderr if stream is None else stream
    if colour is None:
        colour = _colour_supported(stream)
    console = logging.StreamHandler(stream)
    console.setLevel(level)
    console.setFormatter(ConsoleFormatter(colour=colour))
    logger.addHandler(console)

    if log_file is not None:
        path = Path(log_file)
        path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
        file_handler.setLevel(TRACE)
        file_handler.setFormatter(FileFormatter())
        logger.addHandler(file_handler)
        logger.setLevel(min(logger.level, TRACE))

    return logger


# --------------------------------------------------------------------------
# presentation helpers
# --------------------------------------------------------------------------

BANNER = r"""
                ____  _                ____           _
    _ __  _   _|  _ \(_)_   _____ _ __| __ )  ___  __| |
   | '_ \| | | | |_) | \ \ / / _ \ '__|  _ \ / _ \/ _` |
   | |_) | |_| |  _ <| |\ V /  __/ |  | |_) |  __/ (_| |
   | .__/ \__, |_| \_\_| \_/ \___|_|  |____/ \___|\__,_|
   |_|    |___/
"""


def log_banner(logger: logging.Logger | None = None) -> None:
    """Log the pyRiverBed banner, version and citation."""
    logger = logger or get_logger()
    logger.info(BANNER.rstrip("\n"))
    logger.info("   Generate Synthetic Riverbed Topography for Meandering Rivers")
    logger.info("   version %s  |  MIT License", __version__)
    logger.info("   Zhi Li  |  zhil2[at]illinois[dot]edu")
    logger.info("   Li & Garcia (2021), Computers & Geosciences 152, 104755")
    logger.info("")


def format_table(
    rows: Sequence[Sequence[object]],
    headers: Sequence[str] | None = None,
    align: str | Sequence[str] = "left",
    title: str | None = None,
) -> str:
    """Render *rows* as a box-drawing table.

    A tiny replacement for ``tabulate``, which pyRiverBed used to depend on.
    Keeping it in-house is one less package to install.

    Parameters
    ----------
    rows
        The body of the table. Every cell is passed through :func:`str`.
    headers
        Column headers. ``None`` omits the header row.
    align
        ``'left'``, ``'right'`` or ``'center'``, either one value for every
        column or one value per column.
    title
        Optional caption placed above the table.
    """
    body = [[str(cell) for cell in row] for row in rows]
    head = [str(cell) for cell in headers] if headers else None
    ncol = max([len(r) for r in body] + [len(head) if head else 0] or [0])
    if ncol == 0:
        return ""
    for row in body:
        row.extend([""] * (ncol - len(row)))
    if head:
        head.extend([""] * (ncol - len(head)))

    if isinstance(align, str):
        aligns = [align] * ncol
    else:
        aligns = list(align) + ["left"] * (ncol - len(align))

    widths = [
        max([len(row[j]) for row in body] + ([len(head[j])] if head else [0]))
        for j in range(ncol)
    ]

    def pad(text: str, width: int, how: str) -> str:
        if how == "right":
            return text.rjust(width)
        if how == "center":
            return text.center(width)
        return text.ljust(width)

    rule = "+" + "+".join("-" * (w + 2) for w in widths) + "+"
    lines = []
    if title:
        lines.append(title)
    lines.append(rule)
    if head:
        lines.append(
            "| " + " | ".join(pad(head[j], widths[j], "center")
                              for j in range(ncol)) + " |"
        )
        lines.append(rule.replace("-", "="))
    for row in body:
        lines.append(
            "| " + " | ".join(pad(row[j], widths[j], aligns[j])
                              for j in range(ncol)) + " |"
        )
    lines.append(rule)
    return "\n".join(lines)


def log_table(
    rows: Sequence[Sequence[object]],
    headers: Sequence[str] | None = None,
    align: str | Sequence[str] = "left",
    title: str | None = None,
    logger: logging.Logger | None = None,
    level: int = logging.INFO,
) -> None:
    """Log a table built by :func:`format_table`."""
    logger = logger or get_logger()
    logger.log(level, "%s", format_table(rows, headers, align, title))


def format_duration(seconds: float) -> str:
    """Format *seconds* the way a person would say it."""
    if seconds < 1:
        return f"{seconds * 1000:.0f} ms"
    if seconds < 60:
        return f"{seconds:.1f} s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)} min {seconds:.0f} s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)} h {int(minutes)} min"


@contextmanager
def stage(description: str, logger: logging.Logger | None = None,
          level: int = logging.INFO):
    """Announce a stage of work and report how long it took.

    On success the elapsed time is appended; on failure the stage is marked
    as failed and the exception propagates.

    Examples
    --------
    >>> with stage('Calculating curvature'):      # doctest: +SKIP
    ...     curvature = compute()
    """
    logger = logger or get_logger()
    logger.log(level, "+> %s ...", description)
    started = time.perf_counter()
    try:
        yield
    except Exception:
        logger.error("   %s failed after %s", description,
                     format_duration(time.perf_counter() - started))
        raise
    logger.log(level, "   %s done in %s", description,
               format_duration(time.perf_counter() - started))


class ProgressBar:
    """A single-line progress bar for the migration loop.

    Writes directly to the stream with a carriage return rather than going
    through :mod:`logging`, because a bar that scrolls is worse than no bar.
    It disables itself when the stream is not a terminal, so redirected
    output and log files stay clean.
    """

    def __init__(
        self,
        total: int,
        label: str = "",
        stream=None,
        enabled: bool | None = None,
        min_interval: float = 0.1,
    ) -> None:
        self.total = max(int(total), 1)
        self.label = label
        self.stream = sys.stderr if stream is None else stream
        if enabled is None:
            enabled = (
                hasattr(self.stream, "isatty")
                and self.stream.isatty()
                and not os.environ.get("PYRIVERBED_NO_PROGRESS")
            )
        self.enabled = bool(enabled)
        self.min_interval = min_interval
        self._started = time.perf_counter()
        self._last_drawn = 0.0
        self._finished = False

    def update(self, done: int, suffix: str = "") -> None:
        """Redraw the bar for *done* completed items."""
        if not self.enabled or self._finished:
            return
        now = time.perf_counter()
        if done < self.total and now - self._last_drawn < self.min_interval:
            return
        self._last_drawn = now
        fraction = min(max(done / self.total, 0.0), 1.0)
        columns = shutil.get_terminal_size((80, 24)).columns
        elapsed = now - self._started
        eta = elapsed / fraction - elapsed if fraction > 0 else 0.0
        tail = f" {done}/{self.total}"
        if fraction > 0:
            tail += f"  eta {format_duration(eta)}"
        if suffix:
            tail += f"  {suffix}"
        width = max(columns - len(self.label) - len(tail) - 10, 10)
        filled = int(round(width * fraction))
        bar = "#" * filled + "." * (width - filled)
        self.stream.write(
            f"\r{self.label}[{bar}] {fraction * 100:3.0f}%{tail}\033[K"
        )
        self.stream.flush()

    def close(self, message: str | None = None) -> None:
        """Finish the bar, leaving the cursor on a fresh line."""
        if not self.enabled or self._finished:
            self._finished = True
            return
        self._finished = True
        self.stream.write("\r\033[K")
        if message:
            self.stream.write(message + "\n")
        self.stream.flush()

    def __enter__(self) -> "ProgressBar":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()


def log_key_values(
    items: Iterable[tuple[str, object]],
    logger: logging.Logger | None = None,
    level: int = logging.INFO,
    indent: str = "   ",
) -> None:
    """Log ``key: value`` pairs with the keys aligned."""
    logger = logger or get_logger()
    pairs = [(str(k), str(v)) for k, v in items]
    if not pairs:
        return
    width = max(len(k) for k, _ in pairs)
    for key, value in pairs:
        logger.log(level, "%s%s  %s", indent, key.ljust(width), value)
