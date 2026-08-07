"""pyRiverBed -- synthetic riverbed topography for meandering rivers.

Combines the Kinoshita curve with the Beck equations to generate the
equilibrium bed topography of a constant-width meandering river, and steps the
planform forward in time with a linearised bend-theory migration model that
produces both neck and chute cutoffs.

Reference
---------
Li, Z., & Garcia, M. H. (2021). pyRiverBed: A Python framework to generate
synthetic riverbed topography for constant-width meandering rivers.
*Computers & Geosciences*, 152, 104755. doi:10.1016/j.cageo.2021.104755

Quick start
-----------
>>> import pyriverbed as prb                       # doctest: +SKIP
>>> config = prb.default_config()                  # doctest: +SKIP
>>> config.channel.width = 0.8                     # doctest: +SKIP
>>> result = prb.run(config)                       # doctest: +SKIP

Three frontends share this API: the CLI (``pyriverbed``), the GUI
(``pyriverbed gui``) and the notebook helpers in :mod:`pyriverbed.notebook`.
"""

from . import art
from ._version import __version__
from .bed import BedTopography, compute_bed
from .art import ArtStyle
from .config import (ChannelConfig, ChuteCutoffConfig, Config, ConfigError,
                     CurvatureConfig, FlipConfig, KinoshitaConfig, LagConfig,
                     MigrationConfig, NeckCutoffConfig, OutputConfig,
                     default_config, dump_ini, load_config, write_config)
from .logging_utils import setup_logging
from .migration import Cutoff
from .model import RiverBedModel, RunResult, run
from .planform import Centerline

__all__ = [
    "ArtStyle",
    "BedTopography",
    "Centerline",
    "ChannelConfig",
    "ChuteCutoffConfig",
    "Config",
    "ConfigError",
    "CurvatureConfig",
    "Cutoff",
    "FlipConfig",
    "KinoshitaConfig",
    "LagConfig",
    "MigrationConfig",
    "NeckCutoffConfig",
    "OutputConfig",
    "RiverBedModel",
    "RunResult",
    "__version__",
    "art",
    "compute_bed",
    "default_config",
    "dump_ini",
    "load_config",
    "run",
    "setup_logging",
    "write_config",
]
