# pyRiverBed

**Generate Synthetic Riverbed Topography for Meandering Rivers**

[![version](https://img.shields.io/badge/version-2.0.0-blue)](https://github.com/ZhiLiHydro/pyRiverBed/releases)
[![license](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.cageo.2021.104755-orange)](https://www.doi.org/10.1016/j.cageo.2021.104755)

![intro](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_intro.png)

<p align="center">
  <img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_fisk.png" width="100%">
  <br>
  <em>A meander belt grown from an initially straight channel, drawn in the
  style of Harold Fisk's 1944 maps of the lower Mississippi. Every colour is a
  channel course the model passed through; the cream channel is the river
  today. <a href="#art-prints">More prints below.</a></em>
</p>

---

## Contents

* [Publications](#publications) · [Citation](#citation)
* **[Theory guide →](THEORY_GUIDE.md)** — the physics and numerics behind the model
* [What's new in v2.0.0](#whats-new-in-v200) · [Introduction](#introduction) · [Features](#features)
* [Prerequisites](#prerequisites) · [Installation](#installation) · [Usage](#usage)
* [Examples](#examples) · [Art prints](#art-prints) · [Output files](#output-files)
* [Migrating from v1.x](#migrating-from-v1x) · [Testing](#testing)
* [Gallery: v1.x figures and animations](#gallery-figures-and-animations-from-v1x)

---

## Publications

Li, Z., & Garcia, M. H. (2021). pyRiverBed: A Python framework to generate
synthetic riverbed topography for constant-width meandering rivers. Computers &
Geosciences, 152.
doi:[10.1016/j.cageo.2021.104755](https://www.doi.org/10.1016/j.cageo.2021.104755)

Rowley, T., Konsoer, K., Langendoen, E. J., Li, Z., Ursic, M., & Garcia, M. H.
(2021). Relationship of point bar morphology to channel curvature and planform
evolution. Geomorphology, 375.
doi:[10.1016/j.geomorph.2020.107541](https://www.doi.org/10.1016/j.geomorph.2020.107541)

Sun, Y., Song, X., Li, Z., Xu, H., & Bai, Y. (2025). Analytical simulation of
meander morphology from equilibrium to long-term evolution: Impacts of channel
geometry and vegetation-induced coarsening. International Journal of Sediment
Research. doi:[10.1016/j.ijsrc.2025.02.003](https://doi.org/10.1016/j.ijsrc.2025.02.003)

## Citation

```bibtex
@article{li2021pyriverbed,
  title   = {pyRiverBed: A Python framework to generate synthetic riverbed
             topography for constant-width meandering rivers},
  author  = {Li, Zhi and Garcia, Marcelo H.},
  journal = {Computers \& Geosciences},
  volume  = {152},
  pages   = {104755},
  year    = {2021},
  doi     = {10.1016/j.cageo.2021.104755}
}
```

---

## The theory guide

Every equation the model solves, where it comes from, and where it stops being
valid, is written out in **[`THEORY_GUIDE.md`](THEORY_GUIDE.md)**:

| Section | What it covers |
|:---|:---|
| [The idea in one page](THEORY_GUIDE.md#1-the-idea-in-one-page) | Why bed topography is a local function of curvature |
| [Channel planform](THEORY_GUIDE.md#3-channel-planform) | The Kinoshita curve, resampling, Savitzky–Golay smoothing, end extension |
| [Curvature](THEORY_GUIDE.md#4-curvature) | Three estimators, de-spiking, filtering |
| [Curvature phase lag](THEORY_GUIDE.md#5-curvature-phase-lag) | Why the pool sits downstream of the apex |
| [Bed topography](THEORY_GUIDE.md#6-bed-topography) | The scour factor, the Beck profile, and the centerline depth derived from area conservation |
| [1D to 2D](THEORY_GUIDE.md#7-from-1d-centerline-to-2d-channel) | Polyline offsetting at constant width |
| [Mesh and BCs](THEORY_GUIDE.md#8-finite-element-mesh-and-boundary-conditions) | The TELEMAC-ready triangulation |
| [Meander migration](THEORY_GUIDE.md#9-meander-migration) | Linearised bend theory, the sign structure that makes meanders grow, the upstream boundary |
| [Cutoffs](THEORY_GUIDE.md#10-cutoffs) | Neck cutoffs geometrically, chute cutoffs stochastically |
| [Assumptions and limits](THEORY_GUIDE.md#11-assumptions-and-limits-of-validity) | What the model does **not** do, stated plainly |
| [Choosing parameters](THEORY_GUIDE.md#12-choosing-parameters) | A typical range and a calibration route for every input |

If you are going to publish something made with this tool, read that file
first — particularly the assumptions section.

---

## What's new in v2.0.0

v2.0 is a full rewrite of the software around the same science.

**Three frontends, one model**

| | |
|:---|:---|
| **CLI** | `pyriverbed run my_river.ini` — scriptable, with `init`, `convert`, `show` and `art` subcommands |
| **GUI** | `pyriverbed gui` — rebuilt on `ttk`: tabbed, validated, threaded, with a live log and planform preview |
| **Notebook** | `from pyriverbed.notebook import quick_run` — see [`notebooks/pyRiverBed_demo.ipynb`](notebooks/pyRiverBed_demo.ipynb) |

**A readable input file.** The v1.x steering file was a bare column of 39
numbers whose meaning was purely positional. v2 uses a commented, sectioned
key/value file, and every value is validated with a message that names the
parameter. Old steering files still load, and `pyriverbed convert` migrates
them.

```ini
[chute_cutoff]
enabled           = yes    # model chute cutoffs (needs migration)
frequency         = 0.1    # probability per time step; 0.1 = 10% chance each step
entrance          = apex   # apex | inflection
min_sinuosity     = 1.2    # the chute's slope advantage over the bend it replaces
```

**Object-oriented, no globals.** `RiverBedModel` owns its `Config`, so several
models can run in one process — impossible in v1.x, where every parameter was a
module global that Numba baked into its kernels at first compile.

**Fewer dependencies, faster.** Numba and tabulate are gone: the numerics are
vectorised NumPy and the tables are formatted in-house. Only **NumPy, SciPy and
Matplotlib** are required. Neck cutoff detection went from an O(n²) scan of
every node pair to a k-d tree query, which is what makes long runs practical.
Every rewritten kernel was checked against v1.1.0 and agrees to
floating-point round-off.

**Real logging.** Everything goes through `logging`, so a run can be followed
on screen, captured to a file, or piped into the GUI's log pane. Runs report
stage timings, a progress bar with an ETA, and a closing summary.

**Age-shaded meander belts.** Oxbow lakes are now filled according to *when*
they were abandoned, so a belt reads as a stratigraphy instead of an
undifferentiated blob.

**Art prints.** `pyriverbed art` renders a run as a print rather than a figure —
seven styles, led by a homage to Harold Fisk's 1944 maps of the Mississippi
meander belt. [See the gallery.](#art-prints)

### Chute cutoff modeling (from v1.1.0, extended in v2.0)

Besides neck cutoffs, which are geometric and deterministic, pyRiverBed models
**chute cutoffs**, where flow carves a new, shorter channel across the
floodplain while the bend limbs are still well separated. Chute initiation
depends on flood history, bank strength and bar topography — none of which the
model resolves — so it is treated as *conditionally random*: the geometry
decides where a chute is **possible**, and a random draw decides whether one
**happens**. All of it is user-controlled:

| Parameter | Meaning |
|:---|:---|
| `enabled` | Switch chute cutoff modeling on (needs migration to be on) |
| `frequency` | Probability of triggering a chute cutoff in one time step, e.g. `0.1` = 10%. Calibrate as `dt / recurrence_interval` |
| `start_step` | Time steps of spin-up before chute cutoffs are allowed |
| `entrance` | Whether chute channels start and end at bend `apex`es or at `inflection` points |
| `span` | Number of entrance points a chute spans; `2` bypasses one full meander loop |
| `max_valley_angle` | Maximum angle between the chute chord and the valley axis; a chute aligned down-valley takes the largest slope advantage |
| `min_length_widths` | Minimum length of the bypassed reach, in channel widths |
| `min_sinuosity` | Minimum sinuosity of the bypassed reach — the slope advantage itself, and the reason a chute forms at all |
| `end_margin` | Entrance points kept clear of the two ends of the centerline |

---

## Introduction

### For the general public

Meandering rivers erode their outer banks and deposit sediment on their inner
banks. This process makes [point bars](https://en.wikipedia.org/wiki/Point_bar),
which are exposed whenever the water level is relatively low and so are visible
in satellite and aerial imagery. It is interesting and educational to apply this
tool to manufacture your own meandering river at any scale, or to investigate
meandering rivers near you.

### For people working on fluvial geomorphology and earth surface dynamics

This tool can:

* help people working on hydrodynamic and morphodynamic modeling of fluvial
  processes prepare their FEM triangle meshes and boundary condition files;
* help people working on field surveying interpolate bathymetry data in
  unexplored zones during their campaigns;
* help people working on laboratory experiments design their flumes.

## Features

* Two modes:
  * synthetic meandering rivers via the built-in Kinoshita curve calculator;
  * your own real river centerlines, read from file.
* Expands the 1D centerline to a 2D river channel by polyline offsetting.
* Calculates riverbed topography through an analytical method.
* Simulates meander channel migration with a linear bend-theory model.
* Detects **neck cutoffs** geometrically and models **chute cutoffs**
  stochastically, recording every oxbow lake.
* Writes finite element mesh and boundary condition files for
  [TELEMAC](http://www.opentelemac.org/) modeling.
* Renders **art prints** of a run: seven styles, led by a Harold Fisk homage.
* Cross-platform: macOS, Linux and Windows.
* Three interchangeable frontends: CLI, GUI and Jupyter notebook.

## Prerequisites

* Python >= 3.10
* NumPy
* SciPy
* Matplotlib
* tkinter — for the GUI only; it ships with most Python installations
* [imageio](https://imageio.github.io/) — optional, only to write animated GIFs
* [ipywidgets](https://ipywidgets.readthedocs.io/) — optional, only for the
  notebook's interactive form

Using [conda](https://docs.conda.io/) and a virtual environment is recommended:

```bash
conda create -n pyriverbed python numpy scipy matplotlib imageio
conda activate pyriverbed
```

**OR**, with `pip`:

```bash
pip3 install numpy scipy matplotlib imageio
```

On some versions of
[Ubuntu in Windows Subsystem for Linux](https://ubuntu.com/wsl), `tkinter` is
not installed with Python 3; `sudo apt install python3-tk`, or the equivalent on
other platforms, fixes it.

The following two are recommended, but not required:

* [PyInstaller](https://pypi.org/project/PyInstaller/), to freeze the code and
  its dependencies into a single package, i.e. to make executables;
* [Gifsicle](https://www.lcdf.org/gifsicle/), a command-line tool to optimize
  GIFs, and its Python wrapper
  [pygifsicle](https://pypi.org/project/pygifsicle/):

```python
from pygifsicle import optimize
optimize("path_to_my_gif.gif")
```

## Installation

```bash
git clone https://github.com/ZhiLiHydro/pyRiverBed.git
cd pyRiverBed
pip install -e .
```

That puts a `pyriverbed` command on your path. Without installing, everything
still works through `python3 -m pyriverbed` from the repository root.

## Usage

### Command line

```bash
pyriverbed init                  # write a commented pyriverbed.ini to start from
pyriverbed run                   # run it
pyriverbed run my_river.ini -o out   # run a named file, output into out/
pyriverbed convert steering.txt  # migrate a v1.x steering file to the v2 format
pyriverbed show my_river.ini     # echo the configuration with defaults filled in
pyriverbed art my_river.ini      # run, then render art prints of the result
pyriverbed gui                   # open the graphical interface
```

Useful flags: `-n/--steps` to override the number of time steps, `--seed` to
make a stochastic run reproducible, `--dry-run` to validate a configuration
without computing anything, and `-v`/`-q` for more or less output.

A run narrates itself:

```
+> Building channel planform ...
   Building channel planform done in 3 ms
   nodes           1069
   channel length  32.03 m
   sinuosity       3.1632
+> Writing riverbed, bankline and mesh files ...
   wrote kinoshita_topo.xyz  (22449 points)
   wrote kinoshita_mesh.t3s, kinoshita_mesh.dat, ... (22449 nodes, 42720 triangles)
   migrating [########################......]  80%  8000/10000  eta 24 s  sinuosity 2.914, 7 cutoffs
```

### Graphical interface

`pyriverbed gui` opens a tabbed window. Parameters are grouped and validated as
you type, the planform can be previewed before committing to a run, the model
runs on a worker thread so the window stays responsive and can be stopped, and
the log appears in the window rather than in a terminal behind it.

**Planform tab**, previewing a real centerline before the run starts:

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_gui_planform.png" width="100%">

**Cutoffs tab**, with the neck and chute cutoff criteria side by side:

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_gui_cutoffs.png" width="100%">

### Jupyter notebook

[`notebooks/pyRiverBed_demo.ipynb`](notebooks/pyRiverBed_demo.ipynb) walks
through the whole model, from a one-line run to the effect of every parameter:

```python
from pyriverbed.notebook import quick_run, show, cross_section

result = quick_run(n_bends=3, width=0.6, depth=0.15)
show(result)
result.summary()
```

The results are plain objects — `result.bed.z`, `result.cloud`,
`result.centerline.curvature` — so you can take the synthetic bed and do your
own analysis with it.

### Python API

```python
import pyriverbed as prb

config = prb.default_config()
config.channel.width = 0.8
config.migration.enabled = True
config.chute_cutoff.enabled = True
config.chute_cutoff.frequency = 0.05

result = prb.RiverBedModel(config).run()
print(result.summary())
```

## Examples

### Mode 1: generate a Kinoshita curve

Task: reproduce the flume studied in
[this paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008WR007017).

The defaults are exactly this case, so `pyriverbed init` followed by
`pyriverbed run` reproduces it with nothing to change. In the GUI, press
**Run pyRiverBed**; in a notebook, `quick_run()`.

![eg1](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_eg1_kinoshita.png)

The bar–pool asymmetry is the heart of the bed model, and it is worth looking at
directly. The flow depth grows **linearly** towards the outer bank, where the
bed is a scour surface, and decays **exponentially** towards the inner bank,
where it is the depositional face of a point bar. The depth on the centerline is
not free: it follows from requiring the cross-sectional area to stay equal to
`width × depth`, so deepening a pool automatically raises the opposite bar.

![cross section](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_cross_section.png)

### Mode 2: read your own river centerline from file

Task: read centerline coordinates from a file.

The centerline of a randomly picked reach (at 7°32'09.9"S 72°31'16.0"W) of a
randomly picked river — the Juruá River in Brazil — was digitised manually (a
proper centerline extraction tool is recommended for real cases) on a
georeferenced TIFF map, and saved to `jurua.txt`. Set mode to `centerline`,
file name to `jurua.txt`, width to 160 m (estimated), depth to 8 m (arbitrary),
slope to 0 (arbitrary) and lag strength to 6 (estimated); keep the defaults for
everything else.

```ini
[run]
mode            = centerline
centerline_file = example_centerlines/jurua.txt

[channel]
width = 160
depth = 8

[lag]
strength = 6

[curvature]
smoothing_level = 50    # raise until the curvature reads as a bend sequence
```

![eg2](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_eg2_jurua.png)

`smoothing_level` is the parameter that matters most for a real centerline.
Curvature is a second derivative, so digitising noise turns into spurious pools.
Raise it until the curvature signal reads as a sequence of bends rather than
noise — and no further, because over-smoothing flattens the bend apexes and
under-predicts the pool depths.

Eight example centerlines ship in
[`example_centerlines/`](example_centerlines): `jurua.txt`,
`mackey_wabash.txt`, `maier_horseshoe_wabash.txt`, `shimanto.txt`,
`straight.txt`, `trinity.txt`, `ucayali.txt` and `white.txt`.

### Modeling meander migration and cutoffs

Switching migration on steps the planform forward in time. Bends grow,
translate downstream, and eventually cut off — by neck when a loop closes on
itself, or by chute when flow takes a shorter path across the floodplain. Every
abandoned loop is recorded as an oxbow lake and shaded by its age.

```ini
[migration]
enabled = yes
n_steps = 20000
dt      = 86400      # s
e0      = 2e-7       # 1/s -- e0 * dt is the displacement per step, in widths
ub0     = 2          # inlet noise; needed to grow bends from a straight channel

[curvature]
migration_smoothing_level = 1    # keep small, or meanders diffuse away

[neck_cutoff]
enabled = yes

[chute_cutoff]
enabled       = yes
frequency     = 0.02             # a 2% chance per step, i.e. one per ~50 steps
start_step    = 500
min_sinuosity = 1.3
```

![meander belt](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_meander_belt.png)

A run also writes a diagnostics figure. Sinuosity climbs while the bends grow
and drops abruptly at every cutoff; that sawtooth is the signature of a meander
belt in dynamic equilibrium. The mean migration rate is the check that
`e0 × dt` is giving a physically sensible rate.

![diagnostics](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_diagnostics.png)

Synthetic riverbed | Meander belt
:-------------------------:|:-------------------------:
<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_migration_bed.gif"> | <img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_migration_belt.gif">

### What chute cutoffs do to a meander belt

Chute cutoff frequency is the one cutoff parameter with no geometric meaning to
fall back on — it stands in for flood history, bank strength and bar topography,
none of which the model resolves — so it has to be calibrated. This is what it
buys you. Three runs, identical in every respect including the random seed,
differing only in `chute_cutoff.frequency`:

| | `chute_cutoff` | Recurrence at `dt` = 1 day | Result over 12,000 steps |
|:---|:---|:---|:---|
| **Neck only** | `enabled = no` | — | 152 neck, 0 chute — sinuosity **2.96** |
| **Infrequent chutes** | `frequency = 0.0005` | one per ~5 years | 146 neck, 5 chute — sinuosity **2.79** |
| **Frequent chutes** | `frequency = 0.01` | one per ~3 months | 93 neck, 92 chute — sinuosity **2.30** |

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_chute_comparison.png" width="100%">

Reproduce it with:

```bash
pyriverbed run examples/05_chute_frequency.ini -o out_chute
```

Two things are worth reading off the panels.

**Chute cutoffs substitute for neck cutoffs rather than adding to them.** Going
from no chutes to one every three months replaces 59 neck cutoffs with 92 chute
cutoffs; the *total* only rises from 152 to 185. A chute takes a bend out before
its limbs can close, so the neck cutoff that bend was heading for never happens.
The cumulative-cutoff panel shows the two curves separating slowly, not by the
factor you would get if chutes were simply extra events.

**What chutes really change is the ceiling on sinuosity.** With neck cutoffs
alone, every bend grows until it closes on itself, so sinuosity repeatedly
ratchets up towards 4 before collapsing. Frequent chutes truncate bends earlier
and hold sinuosity near 2.5, visibly below the other two runs for the whole run.
The migration rate is almost unaffected: chute frequency does not change how fast
banks move, only how long a bend is allowed to survive.

At a realistic recurrence — the infrequent case, five chutes in the 25 years for
which they were allowed — the effect on the reach-scale statistics is small. That is the honest answer for most
rivers, and it is why `frequency` should be calibrated against an observed cutoff
record rather than guessed: set it to `dt / recurrence_interval`, and run several
seeds before drawing conclusions from any one of them.

## Art prints

A migration run produces, almost incidentally, something people hang on walls.

In 1944 the geologist **Harold Fisk** mapped the meander belt of the lower
Mississippi for the US Army Corps of Engineers. To show every course the river
had taken he drew each one in its own flat colour and let them overlap — the
1880 channel in green, 1820 in salmon-pink, 1765 in light blue, and behind them
the prehistoric courses, a palimpsest on cream paper. The sheets were an
appendix to a dry government report. They are now the most reproduced piece of
fluvial geomorphology there is, sold as prints in a hundred shops.

That is exactly the data pyRiverBed generates: a sequence of channel courses,
and an oxbow lake stamped with the step it was abandoned. So
[`pyriverbed/art.py`](pyriverbed/art.py) draws it the way Fisk drew a survey.

```bash
pyriverbed art examples/04_straight_to_meandering.ini
```

```bash
pyriverbed art my_river.ini -s fisk --paper a2 --art-dpi 300 --format pdf
```

Where [`plotting.py`](pyriverbed/plotting.py) draws *figures* — labelled axes,
colour bars, everything needed to check a number — `art.py` draws *prints*: no
axes, no ticks, no legend. The channel is a filled ribbon of its own width, and
the composition is the whole content. Seven styles ship:

### `fisk` — the meander belt as a palimpsest

Every historical course in its own flat colour, oldest at the back, oxbow lakes
in the colour of the era that abandoned them. The live channel is drawn in
*paper* rather than ink — Fisk's "mighty blank" — so the present river reads as
the hole punched through every course that came before it. Older courses carry
a hatch as well as a hue, which is how Fisk kept overlapping inks apart and what
keeps the print legible in greyscale.

![fisk](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_fisk.png)

### `strata` — one ramp, oldest to newest

The same courses separated by lightness instead of hue, so the print reads as a
single object with depth in time rather than a stack of maps. Overprinting
density stands for how long the river stayed put.

![strata](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_strata.png)

### `blueprint` — cyanotype

Pale channel on Prussian blue with drafting furniture. The grid is spaced in
*channel widths*, which is the scale the model actually thinks in, and it
doubles as the texture.

![blueprint](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_blueprint.png)

### `nocturne` — gold on near-black

The channel stroked several times at increasing width and falling alpha, which
fakes a bloom with no image filter and no extra dependency.

![nocturne](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_nocturne.png)

### `minimal` — one ribbon, a lot of air

![minimal](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_minimal.png)

### `bathymetry` and `contour` — the bed itself

The two styles that show what pyRiverBed actually computes: pools against the
outer banks, point bars against the inner ones. `bathymetry` fills the depth
bands; `contour` keeps the lines and drops the fill.

`bathymetry` | `contour`
:-------------------------:|:-------------------------:
<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_bathymetry.png"> | <img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/v2_art_contour.png">

### From Python or a notebook

```python
from pyriverbed import art
from pyriverbed.notebook import quick_run, show_art, art_styles

result = quick_run("examples/04_straight_to_meandering.ini")

show_art(result, "fisk")                     # inline, screen resolution
art.render(result, "fisk", art.ArtStyle(paper="a2", dpi=300),
           paths=["mississippi.pdf"])        # vector, print resolution
art.save_gallery(result, "prints")           # every style at once
print(art_styles())
```

`ArtStyle` controls the paper (`fit`, the default, shapes the canvas to the
reach; or `a4`/`a3`/`a2`/`letter`/`tabloid`/`square`), the resolution, how many
courses to draw, the cartouche text, the scale bar and the paper texture. Output
format follows the file extension, so `.pdf` and `.svg` give vector prints.

## Output files

| File | Format | Use |
|:---|:---|:---|
| `*_topo.xyz` | 3-column ASCII point cloud | Blue Kenue, GIS, any interpolator |
| `*_boundary.i2s` | closed polyline | the banklines, for Blue Kenue |
| `*_mesh.t3s` | Blue Kenue T3 mesh | build a TELEMAC Selafin geometry file |
| `*_mesh.dat` | Tecplot FETRIANGLE | visualisation |
| `*_BC.cli` | TELEMAC boundary conditions | use directly in a TELEMAC run |
| `*_BC.bc2` | Blue Kenue BC | inspect and edit BC codes and metadata |
| `*_pyriverbed.png/pdf` | figure | the three-panel summary |
| `*_diagnostics.png` | figure | sinuosity and migration rate time series |
| `*_sinuosity.txt`, `*_mean_migration_rate.txt` | ASCII | one value per time step |
| `*_migration0.gif`, `*_migration1.gif` | animation | the bed and the meander belt |
| `art/*_<style>.png` | print | art prints, from `pyriverbed art` |
| `pyriverbed.log` | text | the full run log |

## Documentation

* **[`THEORY_GUIDE.md`](THEORY_GUIDE.md)** — the physics and the numerics, every
  governing equation, and an explicit statement of the assumptions and their
  limits. [Section index above.](#the-theory-guide)
* **[`notebooks/pyRiverBed_demo.ipynb`](notebooks/pyRiverBed_demo.ipynb)** — a
  guided tour of every parameter, in 13 sections.
* **[`examples/`](examples)** — five commented input files, and a note on the
  handful of parameters that interact.
* **`pyriverbed.ini`** — written by `pyriverbed init`, with a one-line
  explanation beside every parameter.

## Migrating from v1.x

| v1.x | v2.0 |
|:---|:---|
| `python3 gui4pyriverbed.py` | `pyriverbed gui` |
| `python3 pyriverbed.py` | `pyriverbed run` |
| `steering.txt` (39 positional numbers) | `pyriverbed.ini` (commented sections) |
| — | `pyriverbed convert steering.txt` |

Old steering files are read as they are, so nothing has to be converted to keep
working. A v1.x file is loaded in a compatibility mode that reproduces v1.x
output; when you convert one to the v2 format, three settings are worth bringing
up to the v2 defaults:

* `curvature.migration_smoothing_level = 2`. A v1.x file smooths as many times
  per migration step as it does when building the centerline, which for a long
  run means the smoother flattens the meanders faster than the model grows them.
* `migration.end_taper_widths = 2` and `neck_cutoff.end_margin_widths = 2`.
  These keep the reach ends, where the artificial straight extensions and the
  inlet transient live, out of the migration and out of the cutoff search.
* `flip.transverse = no`. The converter sets this for you; it is now a purely
  cosmetic mirror of the bed about the centerline.

## Testing

```bash
python -m unittest discover -s tests -v
```

84 tests, using only the standard library's `unittest`, covering the
configuration round-trip and legacy reader, every curvature estimator against
an analytical circle, the bed model's area conservation, polyline offsetting
against analytical offsets, mesh connectivity, both cutoff mechanisms, the end
taper, every art style, the CLI and notebook APIs, and the LaTeX in these docs
against the constructs GitHub's Markdown renderer silently breaks.

## License

[MIT License](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/LICENSE)

---

## Gallery: figures and animations from v1.x

These were produced with pyRiverBed v1.x and are kept here as a record of the
project. The physics is the same; only the plotting and the frontends have
changed.

### The v1.x tkinter GUI

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_gui1.gif">

### Mode 1 and Mode 2 example figures

Mode 1: Kinoshita curve, annotated with the Blue Kenue and TELEMAC output

![v1 eg1](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_eg1.png)

Mode 2: the Juruá River, with the Landsat imagery it was digitised from

![v1 eg2](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_eg2.png)

### Meander migration animations

Kinoshita curve — synthetic riverbed | Kinoshita curve — river centerline
:-------------------------:|:-------------------------:
<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/kinoshita_migration0.gif"> | <img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/kinoshita_migration1.gif">

Juruá River — synthetic riverbed | Juruá River — river centerline
:-------------------------:|:-------------------------:
<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/jurua_migration0.gif"> | <img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/jurua_migration1.gif">

### Long-term evolution

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/longterm_a.gif">

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/longterm_b.gif">
