"""Command line interface for pyRiverBed.

Subcommands
-----------
``run``      generate topography from an input file (the default)
``init``     write a commented input file to start from
``convert``  turn a v1.x ``steering.txt`` into a v2 input file
``show``     print a configuration as pyRiverBed understands it
``art``      run, then render art prints of the result
``gui``      launch the graphical interface

Run ``pyriverbed --help`` or ``pyriverbed <subcommand> --help`` for details.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from ._version import __version__
from .config import (Config, ConfigError, default_config, dump_ini,
                     load_config, write_config)
from .logging_utils import get_logger, log_banner, setup_logging

log = get_logger(__name__)

DEFAULT_INPUT = "pyriverbed.ini"

_EPILOG = """\
examples:
  pyriverbed init                        write a commented pyriverbed.ini
  pyriverbed run                         run pyriverbed.ini
  pyriverbed run my_river.ini -o out     run and write everything into out/
  pyriverbed run steering.txt            run a v1.x steering file as-is
  pyriverbed convert steering.txt        migrate a v1.x file to the v2 format
  pyriverbed show my_river.ini           echo the configuration, defaults filled
  pyriverbed art my_river.ini            run, then render every art print
  pyriverbed art -s fisk --paper a2      one style, at poster size
  pyriverbed gui                         open the graphical interface

Documentation: https://github.com/ZhiLiHydro/pyRiverBed
Theory:        THEORY_GUIDE.md
"""


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="pyriverbed",
        description="Generate synthetic riverbed topography for meandering "
                    "rivers.",
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version",
                        version=f"pyRiverBed {__version__}")

    common = argparse.ArgumentParser(add_help=False)
    verbosity = common.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="count", default=0,
                           help="more output; repeat for debug detail")
    verbosity.add_argument("-q", "--quiet", action="store_true",
                           help="warnings and errors only")
    common.add_argument("--no-color", action="store_true",
                        help="disable coloured output")
    common.add_argument("--log-file", metavar="PATH",
                        help="write a log file here, overriding the input file")

    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser(
        "run", parents=[common], help="generate topography (default)",
        description="Generate riverbed topography from an input file.")
    run_parser.add_argument("input", nargs="?", default=DEFAULT_INPUT,
                            help=f"input file (default: {DEFAULT_INPUT})")
    run_parser.add_argument("-o", "--output", metavar="DIR",
                            help="output directory, overriding the input file")
    run_parser.add_argument("-n", "--steps", type=int, metavar="N",
                            help="number of migration time steps, overriding "
                                 "the input file")
    run_parser.add_argument("--seed", type=int, metavar="N",
                            help="seed the random generator for a reproducible "
                                 "stochastic run")
    run_parser.add_argument("--dry-run", action="store_true",
                            help="validate and report the configuration "
                                 "without computing anything")

    init_parser = subparsers.add_parser(
        "init", parents=[common], help="write a commented input file",
        description="Write a fully commented input file holding the defaults.")
    init_parser.add_argument("output", nargs="?", default=DEFAULT_INPUT,
                             help=f"file to write (default: {DEFAULT_INPUT})")
    init_parser.add_argument("-f", "--force", action="store_true",
                             help="overwrite an existing file")

    convert_parser = subparsers.add_parser(
        "convert", parents=[common],
        help="convert a v1.x steering file to the v2 format",
        description="Read a v1.x steering.txt and write the equivalent v2 "
                    "input file.")
    convert_parser.add_argument("input", help="v1.x steering file")
    convert_parser.add_argument("output", nargs="?",
                                help=f"file to write "
                                     f"(default: {DEFAULT_INPUT})")
    convert_parser.add_argument("-f", "--force", action="store_true",
                                help="overwrite an existing file")

    show_parser = subparsers.add_parser(
        "show", parents=[common], help="print a configuration",
        description="Print a configuration in the v2 format, with every "
                    "default filled in.")
    show_parser.add_argument("input", nargs="?",
                             help="input file; omit to show the defaults")

    art_parser = subparsers.add_parser(
        "art", parents=[common], help="render art prints from a run",
        description="Run the model and render art prints of the result: "
                    "no axes, no legend, just the river on a flat ground. "
                    "The 'fisk' style is a homage to Harold Fisk's 1944 maps "
                    "of the Mississippi meander belt, and needs migration to "
                    "be on so that there are historical courses to draw.")
    art_parser.add_argument("input", nargs="?", default=DEFAULT_INPUT,
                            help=f"input file (default: {DEFAULT_INPUT})")
    art_parser.add_argument("-o", "--output", metavar="DIR",
                            help="output directory, overriding the input file")
    art_parser.add_argument("-n", "--steps", type=int, metavar="N",
                            help="number of migration time steps, overriding "
                                 "the input file")
    art_parser.add_argument("--seed", type=int, metavar="N",
                            help="seed the random generator for a reproducible "
                                 "stochastic run")
    art_parser.add_argument("-s", "--style", action="append", metavar="NAME",
                            help="style to render; repeat for several. "
                                 "Default: every style")
    art_parser.add_argument("--paper", default="a3", metavar="SIZE",
                            help="paper size: a4, a3, a2, letter, tabloid, "
                                 "square (default: a3)")
    art_parser.add_argument("--art-dpi", type=int, default=300, metavar="N",
                            help="output resolution (default: 300)")
    art_parser.add_argument("--format", action="append", metavar="EXT",
                            help="output format; repeat for several. "
                                 "png, pdf or svg (default: png)")
    art_parser.add_argument("--title", help="cartouche title; '' to omit")
    art_parser.add_argument("--subtitle", help="cartouche subtitle")
    art_parser.add_argument("--no-texture", action="store_true",
                            help="switch off the paper mottling")
    art_parser.add_argument("--list-styles", action="store_true",
                            help="list the styles and exit")

    subparsers.add_parser(
        "gui", parents=[common], help="launch the graphical interface",
        description="Open the pyRiverBed graphical interface.")
    return parser


def _configure_logging(args: argparse.Namespace,
                       log_file: str | Path | None = None) -> None:
    """Set up logging from the parsed verbosity flags."""
    if getattr(args, "quiet", False):
        level = logging.WARNING
    elif getattr(args, "verbose", 0) >= 2:
        level = logging.DEBUG
    elif getattr(args, "verbose", 0) == 1:
        level = 15  # TRACE
    else:
        level = logging.INFO
    setup_logging(level=level, log_file=log_file,
                  colour=False if getattr(args, "no_color", False) else None,
                  stream=sys.stdout)


def _resolve_log_file(args: argparse.Namespace,
                      config: Config | None) -> Path | None:
    """Decide where the log file goes, if anywhere."""
    if getattr(args, "log_file", None):
        return Path(args.log_file)
    if config is not None and config.output.log_file:
        return config.output.path / config.output.log_file
    return None


def command_run(args: argparse.Namespace) -> int:
    """Handle ``pyriverbed run``."""
    _configure_logging(args)
    path = Path(args.input)
    if not path.exists():
        log.error("input file '%s' not found", path)
        log.error("write one with:  pyriverbed init %s", path)
        return 2
    config = load_config(path)

    if args.output:
        config.output.directory = args.output
    if args.steps is not None:
        config.migration.n_steps = args.steps
        config.migration.enabled = args.steps > 0
    if args.seed is not None:
        config.migration.seed = args.seed
    config.validate()

    # Re-open logging now that the output directory is known.
    _configure_logging(args, _resolve_log_file(args, config))
    log_banner()

    if args.dry_run:
        from .model import RiverBedModel
        log.info("dry run: validating configuration only")
        RiverBedModel(config).log_configuration()
        log.info("configuration is valid")
        return 0

    from .model import run as run_model
    result = run_model(config)
    return 0 if result.steps_completed >= 0 else 1


def command_art(args: argparse.Namespace) -> int:
    """Handle ``pyriverbed art``."""
    from . import art as art_module

    _configure_logging(args)
    if args.list_styles:
        for name, function in sorted(art_module.STYLES.items()):
            summary = (function.__doc__ or "").strip().splitlines()[0]
            log.info("  %-11s %s", name, summary)
        return 0

    unknown = [name for name in (args.style or [])
               if name not in art_module.STYLES]
    if unknown:
        log.error("unknown art style(s): %s", ", ".join(unknown))
        log.error("available: %s", ", ".join(sorted(art_module.STYLES)))
        return 2
    if args.paper not in art_module.PAPERS:
        log.error("unknown paper size '%s'; available: %s", args.paper,
                  ", ".join(sorted(art_module.PAPERS)))
        return 2

    path = Path(args.input)
    if not path.exists():
        log.error("input file '%s' not found", path)
        log.error("write one with:  pyriverbed init %s", path)
        return 2
    config = load_config(path)
    if args.output:
        config.output.directory = args.output
    if args.steps is not None:
        config.migration.n_steps = args.steps
        config.migration.enabled = args.steps > 0
    if args.seed is not None:
        config.migration.seed = args.seed
    config.validate()

    _configure_logging(args, _resolve_log_file(args, config))
    log_banner()
    if not config.migration.enabled:
        log.warning("migration is off, so there is only one channel course; "
                    "the 'fisk' and 'strata' styles have nothing to layer")

    from .model import run as run_model
    result = run_model(config)

    style = art_module.ArtStyle(
        paper=args.paper, dpi=args.art_dpi, title=args.title,
        subtitle=args.subtitle, texture=0.0 if args.no_texture else 0.35)
    log.info("")
    log.info("Rendering art prints")
    art_module.save_gallery(result, config.output.path / "art",
                            prefix=result.config.name or "pyriverbed",
                            styles=args.style, style=style,
                            formats=tuple(args.format or ("png",)))
    return 0


def command_init(args: argparse.Namespace) -> int:
    """Handle ``pyriverbed init``."""
    _configure_logging(args)
    path = Path(args.output)
    if path.exists() and not args.force:
        log.error("'%s' already exists; pass --force to overwrite", path)
        return 2
    write_config(default_config(), path)
    log.info("edit it, then run:  pyriverbed run %s", path)
    return 0


def command_convert(args: argparse.Namespace) -> int:
    """Handle ``pyriverbed convert``."""
    _configure_logging(args)
    source = Path(args.input)
    target = Path(args.output or DEFAULT_INPUT)
    if not source.exists():
        log.error("input file '%s' not found", source)
        return 2
    if target.exists() and not args.force:
        log.error("'%s' already exists; pass --force to overwrite", target)
        return 2
    config = load_config(source)
    write_config(config, target)
    log.info("converted %s -> %s", source, target)
    return 0


def command_show(args: argparse.Namespace) -> int:
    """Handle ``pyriverbed show``."""
    _configure_logging(args)
    config = load_config(args.input) if args.input else default_config()
    sys.stdout.write(dump_ini(config))
    return 0


def command_gui(args: argparse.Namespace) -> int:
    """Handle ``pyriverbed gui``."""
    _configure_logging(args)
    try:
        from .gui import main as gui_main
    except ImportError as exc:
        log.error("cannot start the GUI: %s", exc)
        log.error("tkinter is part of the standard library but is packaged "
                  "separately on some Linux distributions; try "
                  "'sudo apt install python3-tk'")
        return 3
    return gui_main()


_COMMANDS = {
    "run": command_run,
    "init": command_init,
    "convert": command_convert,
    "show": command_show,
    "art": command_art,
    "gui": command_gui,
}


def main(argv: list[str] | None = None) -> int:
    """Entry point for the ``pyriverbed`` command.

    Parameters
    ----------
    argv
        Argument list. ``None`` uses :data:`sys.argv`.

    Returns
    -------
    int
        Process exit status.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    # Make 'run' the default subcommand so that both `pyriverbed` and
    # `pyriverbed my_river.ini` do the obvious thing.
    if argv and argv[0] not in _COMMANDS and not argv[0].startswith("-"):
        argv.insert(0, "run")
    elif not argv:
        argv = ["run"]
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0

    try:
        return _COMMANDS[args.command](args)
    except ConfigError as exc:
        log.error("%s", exc)
        return 2
    except FileNotFoundError as exc:
        log.error("%s", exc)
        return 2
    except KeyboardInterrupt:
        log.warning("interrupted")
        return 130


if __name__ == "__main__":
    sys.exit(main())
