# pyRiverBed
Generate Synthetic Riverbed Topography for Meandering Rivers

![intro](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_intro.png)

## What's new in v1.1.0

Chute cutoff modeling, introduced in v1.0.1, is now fully user-controlled from the steering file instead of being hard-coded.

Besides the meander neck cutoff submodel, which detects cutoffs geometrically, chute cutoffs are triggered stochastically: at every time step of a meander migration run, a chute cutoff happens at a user-given **frequency**, and, once triggered, the chute channel is picked randomly among all the chute channels that satisfy the user-given geometric criteria. The chute channel replaces the channel reach it bypasses, and the bypassed reach becomes an oxbow lake. The following chute cutoff parameters are now available in the GUI and in the steering file:

| Parameter | Meaning |
|:---|:---|
| Chute cutoff | Switch chute cutoff modeling on or off (needs meander migration to be on) |
| Frequency | Probability of triggering a chute cutoff in one time step, e.g. 0.1 = 10% |
| Starting step | Number of time steps of spin-up before chute cutoffs are allowed |
| Entrance at | Whether chute channels start and end at bend apexes or at inflection points |
| Span | Number of entrance points a chute channel spans, e.g. 2 = one full meander bend |
| Max angle | Maximum angle, in degrees, between the chute channel and the valley axis; a chute channel well aligned with the down-valley direction takes a larger slope advantage over the flow path along the channel |
| Min length | Minimum length of the bypassed channel reach, normalized by channel width |
| Margin at ends | Number of entrance points kept away from the two ends of the centerline |

Steering files written by earlier versions of the GUI are still readable: the chute cutoff parameters simply fall back to their default values, and chute cutoff modeling stays off.

## Introduction

### For general public: 

Meandering rivers erode their outer banks and deposit sediments on their inner banks. This process makes [point bars](https://en.wikipedia.org/wiki/Point_bar), which are always exposed when water level is relatively low and thus visible in satellite/aerial imageries. It is interesting and educational to apply this tool to manufacture your own meandering river at any scales, or to investigate meandering rivers near you. 

### For folks working on fluvial geomorphology and earth surface dynamic processes:

This tool could:

* help people working on hydrodynamic and morphodynamic modeling of fluvial processes with preparing their FEM triangle meshes and boundary condition files

* help people working on field surveying with interpolating bathymetry data in unexplored zones during their campaigns

* help people working on laboratory experiments with designing their flumes

## Features

* Provides two modes: 
  * Making synthetic meandering rivers via the built-in Kinoshita Curve calculator
  * Reading users' own real river centerlines 

* Expands 1D centerline to 2D river channel by polyline offsetting

* Calculates riverbed topography through an analytical method

* Simulates meander channel migration using linear channel migration model

* Makes finite element mesh files and boundary condition files for [TELEMAC](http://www.opentelemac.org/) numerical modeling

* Cross-platform, runs on macOS, Linux, Windows machines

* Uses a tkinter-based graphical user interface (GUI) as front end

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_gui1.gif">

## Prerequisites

* Python >= 3.5
* tkinter (in most cases it comes with Python installation)
* Numpy
* Scipy
* Matplotlib
* [Numba](http://numba.pydata.org/)
* [tabulate](https://pypi.org/project/tabulate/)
* [imageio](https://imageio.github.io/) (only needed if simulating meander channel migration)

Using [Anaconda](https://www.anaconda.com/distribution/) and using `conda` to create a virtual environment are recommended:

```
conda create -n pyriverbed python numpy scipy matplotlib numba tabulate imageio
conda activate pyriverbed
```

**OR**, if `pip` is more preferred, just install `Numpy`, `Scipy`, `Matplotlib`, `Numba`, `tabulate` and `imageio` by:

```
pip3 install numpy scipy matplotlib numba tabulate imageio
```

In some versions of [Ubuntu in Windows Subsystem for Linux](https://ubuntu.com/wsl), `tkinter` may not be installed with Python3, then installation of tkinter by `sudo apt install python3-tk` in Ubuntu or a similar command in other platforms is required.

The following two are recommended, but not required:

* [PyInstaller](https://pypi.org/project/PyInstaller/) (to freeze python codes and dependencies into a single package, i.e., making executables);
* [Gifsicle](https://www.lcdf.org/gifsicle/) (a command-line tool to optimize GIFs) and its Python wrapper [pygifsicle](https://pypi.org/project/pygifsicle/). 

Example:

```
from pygifsicle import optimize
optimize("path_to_my_gif.gif")
```

## Installation and usage

### Step 0: Download codes

`git clone https://github.com/ZhiLiHydro/pyRiverBed.git`

### Step 1: Start the GUI to prepare the steering file

`python3 gui4pyriverbed.py` first, then enter model parameters, and click `Generate steering file` button to prepare the steering file `steering.txt`

### Step 2: Run pyRiverBed

GUI method: click `Run pyRiverBed` button in GUI to run pyRiverBed

CLI method: quit the GUI and then type `python3 pyriverbed.py` in command line to run pyRiverBed

Both methods to run pyRiverBed work with zero differences

## Examples

### Mode 1: Generate Kinoshita Curve

Task: Reproduce the flume studied in [this paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008WR007017).

The default parameters in GUI are pre-typed for this case, so nothing needs to be changed. Click `Generate steering file` button, then click `Run pyRiverBed` button to run and check the results. 

![intro](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_eg1.png)

#### Modeling meander migration

Synthetic riverbed | River centerline
:-------------------------:|:-------------------------:
<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/kinoshita_migration0.gif">  |  <img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/kinoshita_migration1.gif">

### Mode 2: Read your own river centerline from file

Task: Read centerline coordinates from file. 

The river centerline of a randomly picked reach (at 7°32'09.9"S 72°31'16.0"W) of a randomly picked river (it's Juruá River in Brazil) is discretized manually (well, a fancy centerline extraction tool is recommended here for real cases) on a georeferenced TIFF map. The river centerline coordinates is saved to `jurua.txt`. Open the GUI to type in the followings: mode = 2, file name = 'jurua.txt', width = 160 m (estimated), depth = 8 m (arbitrary), slope = 0 (arbitrary), lag strength = 6 (estimated), flip in transverse direction = no. Keep the defaults for other parameters. Click `Generate steering file` button, then click `Run pyRiverBed` button in GUI to run and see the results.

![intro](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/pyRiverBed_eg2.png)

#### Modeling meander migration

Synthetic riverbed | River centerline
:-------------------------:|:-------------------------:
<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/jurua_migration0.gif">  |  <img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/jurua_migration1.gif">

## Related publications

Li, Z., & Garcia, M. H. (2021). pyRiverBed: A Python framework to generate synthetic riverbed topography for constant-width meandering rivers. Computers & Geosciences, 152. doi:[10.1016/j.cageo.2021.104755](https://www.doi.org/10.1016/j.cageo.2021.104755)

Rowley, T., Konsoer, K., Langendoen, E. J., Li, Z., Ursic, M., & Garcia, M. H. (2021). Relationship of point bar morphology to channel curvature and planform evolution. Geomorphology, 375. doi:[10.1016/j.geomorph.2020.107541](https://www.doi.org/10.1016/j.geomorph.2020.107541)

Sun, Y., Song, X., Li, Z., Xu, H., & Bai, Y. (2025). Analytical simulation of meander morphology from equilibrium to long-term evolution: Impacts of channel geometry and vegetation-induced coarsening. International Journal of Sediment Research. doi:[10.1016/j.ijsrc.2025.02.003](https://doi.org/10.1016/j.ijsrc.2025.02.003)

## License

[MIT License](https://github.com/ZhiLiHydro/pyRiverBed/blob/master/LICENSE)

## More GIFs

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/longterm_a.gif">

<img src="https://github.com/ZhiLiHydro/pyRiverBed/blob/master/img/longterm_b.gif">
