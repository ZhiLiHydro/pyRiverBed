# pyRiverBed Theory Guide

This document describes the physics and the numerics behind pyRiverBed. It is
the companion of the [README](README.md): the README explains *how to run* the
model, this guide explains *what the model computes and why*.

The reference for the framework is:

> Li, Z., & Garcia, M. H. (2021). pyRiverBed: A Python framework to generate
> synthetic riverbed topography for constant-width meandering rivers.
> *Computers & Geosciences*, 152, 104755.
> doi:[10.1016/j.cageo.2021.104755](https://www.doi.org/10.1016/j.cageo.2021.104755)

Everything below is implemented in the `pyriverbed` package. Where a section
maps onto a specific module, the module is named in the heading.

---

## Contents

1. [The idea in one page](#1-the-idea-in-one-page)
2. [Notation](#2-notation)
3. [Channel planform](#3-channel-planform)
4. [Curvature](#4-curvature)
5. [Curvature phase lag](#5-curvature-phase-lag)
6. [Bed topography](#6-bed-topography)
7. [From 1D centerline to 2D channel](#7-from-1d-centerline-to-2d-channel)
8. [Finite element mesh and boundary conditions](#8-finite-element-mesh-and-boundary-conditions)
9. [Meander migration](#9-meander-migration)
10. [Cutoffs](#10-cutoffs)
11. [Assumptions and limits of validity](#11-assumptions-and-limits-of-validity)
12. [Choosing parameters](#12-choosing-parameters)
13. [References](#13-references)

---

## 1. The idea in one page

A meandering river with a constant width has, to a very good first
approximation, a bed topography that is a *local function of the channel
curvature*. Flow around a bend is superelevated on the outer bank; the
resulting cross-stream pressure gradient is not balanced by the cross-stream
component of the near-bed velocity, so a secondary (helical) flow cell
develops. It sweeps sediment from the outer bank towards the inner bank,
scouring a **pool** against the outer bank and building a **point bar**
against the inner bank. The steeper the bend, the stronger the secondary flow
and the steeper the resulting transverse bed slope.

pyRiverBed turns that statement into a generator:

```
        planform                curvature              bed topography
   ┌────────────────┐      ┌────────────────┐      ┌────────────────┐
   │ Kinoshita curve│      │  C(s) = dθ/ds  │      │  transverse    │
   │       or       │─────>│  filtered      │─────>│  slope S_T     │─────> z(s,n)
   │ your own x,y   │      │  phase-lagged  │      │  Beck profile  │
   └────────────────┘      └────────────────┘      └────────────────┘
                                   │
                                   │  linear bend theory
                                   v
                          ┌────────────────┐
                          │ near-bank      │
                          │ excess velocity│─────> bank migration ─────> cutoffs
                          └────────────────┘
```

Three things are worth emphasising because they are what makes the result look
like a real river rather than a textbook sketch:

* **The bed is not symmetric.** On the pool side the flow depth grows
  *linearly* away from the centerline; on the bar side it decays
  *exponentially*. The two are stitched together by a mass-conservation
  constraint (§6).
* **The bed lags the curvature.** The secondary flow needs a finite distance
  to spin up and to spin down, so the deepest point of a bend sits
  *downstream* of the bend apex. pyRiverBed reproduces this with an
  upstream-weighted moving average of the curvature signal (§5).
* **The planform can evolve.** The same curvature field drives a linearised
  bend-theory migration model, which grows meanders, translates them
  downstream, and eventually cuts them off (§9, §10).

---

## 2. Notation

| Symbol | Meaning | Unit | Input name |
|:---|:---|:---|:---|
| $s$ | streamwise (curvilinear) coordinate along the centerline | m | |
| $n$ | transverse coordinate, $n \in [-b, b]$, positive to the left of the flow | m | |
| $L$ | total centerline length | m | |
| $B$ | channel width (constant) | m | `width` |
| $b = B/2$ | channel half-width | m | |
| $H$ | reach-averaged flow depth | m | `depth` |
| $\beta = b/H$ | half-width-to-depth ratio | – | |
| $I$ | longitudinal (channel) bed slope | – | `slope` |
| $C(s)$ | centerline curvature, positive for a left turn | m⁻¹ | |
| $\theta(s)$ | angle between the local tangent and the $x$ axis | rad | |
| $\lambda$ | arc wavelength of the Kinoshita curve | m | `arc_wavelength` |
| $\theta_0$ | maximum angular amplitude | rad | `max_angular_amplitude` |
| $J_s, J_f$ | Kinoshita skewness and flatness coefficients | – | `skewness`, `flatness` |
| $\Delta s$ | streamwise node spacing | m | `ds` |
| $N$ | number of polyline offsets per side | – | `n_offsets` |
| $\Delta n = b/N$ | transverse node spacing | m | |
| $A$ | scour factor | – | |
| $S_T$ | transverse bed slope | – | |
| $\xi_{S_T}$ | transverse slope corrector | – | `transverse_slope_corrector` |
| $h_c$ | flow depth on the centerline | m | |
| $\zeta(s,n)$ | flow depth | m | |
| $z(s,n)$ | bed elevation | m | |
| $\ell$ | phase lag strength, in channel widths | – | `lag_strength` |
| $u_b$ | dimensionless near-bank excess velocity | – | |
| $C_{f0}$ | reach-averaged friction coefficient | – | `cf0` |
| $Fr_0$ | reach-averaged Froude number | – | `fr0` |
| $E_0$ | bank erosion coefficient | s⁻¹ | `e0` |
| $\Delta t$ | migration time step | s | `dt` |
| $\Omega$ | sinuosity, $L / \Vert P_L - P_0 \Vert$ | – | |

---

## 3. Channel planform

*Module: `pyriverbed.planform`*

### 3.1 Mode 1 — the Kinoshita curve

The Kinoshita curve (Kinoshita, 1961; Parker et al., 1982; Abad & Garcia,
2009) is the standard idealised meander planform. It is defined not by its
coordinates but by its **direction angle** as a function of arc length:

```math
\theta(s) = \theta_0 \sin\left(\frac{2\pi s}{\lambda}\right)
 + \theta_0^{3}\left[ J_s \cos\left(\frac{6\pi s}{\lambda}\right)
 - J_f \sin\left(\frac{6\pi s}{\lambda}\right)\right]
```

The first harmonic is a **sine-generated curve** (Langbein & Leopold, 1966) —
the shape a river takes if it minimises the variance of its direction changes.
The third-harmonic correction is what makes the curve interesting:

* $J_s$ (**skewness**) makes the bends asymmetric in the streamwise direction,
  producing the upstream-skewed, "fattened-downstream" bends that real rivers
  show. It is the term that breaks the up-/downstream symmetry.
* $J_f$ (**flatness**) flattens the bend apexes and sharpens the crossings,
  producing the boxy, high-amplitude bends of tightly meandering rivers.

With $J_s = J_f = 0$ the curve degenerates to the sine-generated curve.

Coordinates follow by integrating the tangent:

```math
x(s) = \int_0^s \cos\theta \mathrm{d}s', \qquad
y(s) = \int_0^s \sin\theta \mathrm{d}s'
```

which pyRiverBed evaluates as a running rectangle sum on a uniform grid of
spacing $\Delta s$. The curve is generated over `n_bends` arc wavelengths.

Because $\theta$ is prescribed analytically, the curvature is available
analytically too:

```math
C(s) = \frac{\mathrm{d}\theta}{\mathrm{d}s}
```

evaluated by first differences on the same grid.

### 3.2 Mode 2 — your own centerline

In Mode 2 the model reads an $(x, y)$ polyline from a text file: two columns,
one node per line, in any projected (metric) coordinate system. Digitised or
extracted centerlines are almost never uniformly spaced, and both the
curvature estimate and the offsetting step assume near-uniform spacing, so the
polyline is immediately **resampled** (§3.3) and **smoothed** (§3.4).

### 3.3 Resampling

The arc length of the input polyline is computed cumulatively,

```math
s_0 = 0, \qquad s_j = s_{j-1} + \sqrt{(x_j - x_{j-1})^2 + (y_j - y_{j-1})^2}
```

and the polyline is re-interpolated linearly onto
$N_r + 1$ equally spaced stations with

```math
N_r = \left\lfloor \frac{s_{\text{end}}}{\Delta n} \right\rfloor + 1,
\qquad \Delta n = \frac{b}{N}
```

i.e. **the streamwise node spacing is set equal to the transverse node
spacing**. This is deliberate: it makes the resulting point cloud
quasi-isotropic, which is what you want if the output is going to become a
finite element mesh (§8) — triangles come out close to equilateral rather than
as slivers.

The two end nodes are pinned so that resampling never shortens the reach.

### 3.4 Smoothing

Curvature is a second derivative of position, so it amplifies noise
brutally: a digitising error of a fraction of a channel width produces a
curvature spike that would put a spurious pool in the middle of a straight
reach. pyRiverBed therefore passes the coordinates through a
**Savitzky–Golay filter** — a 5-point, second-order polynomial least-squares
filter, applied to $x$ and $y$ independently, repeatedly.

The number of passes is set by the smoothing level $m$:

```math
n_{\text{passes}} =
\begin{cases}
m & m < 39\\[2pt]
\mathrm{round}\left(1.1^{ m}\right) & m \geq 39
\end{cases}
```

The exponential branch exists purely to give the GUI slider a useful dynamic
range: levels 0–38 are "gentle", and levels above 38 escalate quickly to the
very heavy smoothing that noisy satellite-derived centerlines sometimes need.
The two end nodes are restored after every pass so that smoothing cannot
retract the reach.

A Savitzky–Golay filter is the right choice here (rather than a moving
average) because it preserves the *amplitude and position of extrema* — the
bend apexes — while removing high-wavenumber noise. A moving average would
flatten the apexes and therefore systematically under-predict pool depth.

**Two smoothing levels, for two different jobs.** Smoothing is applied in two
places, and they need very different strengths:

* `smoothing_level` is applied **once**, when the centerline is built. Its job
  is to clean up digitising noise in an imported centerline, and it can be
  large — tens or hundreds of passes for a noisy satellite-derived line.
* `migration_smoothing_level` is applied **after every migration time step**.
  Its job is only to remove the small-scale roughness that pointwise bank
  displacement introduces, and it must be **small**, one or two passes.

The reason is that smoothing is diffusive. Even a Savitzky–Golay filter, applied
thousands of times, erodes finite-amplitude bends. If the diffusion it produces
per step exceeds the growth the migration model produces per step, the reach
does not meander at all: it decays monotonically towards a straight line. The
two rates are independent, so the only way to keep the model in the growing
regime is to keep the per-step smoothing much weaker than the once-off
smoothing.

> **Changed in v2.0.** v1.x had a single level and used it in both places, so a
> migration run started from the GUI defaults (20 passes) flattened itself. A
> v1.x steering file still reproduces v1.x behaviour, with a warning.

### 3.5 End extension

Both the offsetting algorithm (§7) and the migration model (§9) behave badly
at the ends of the reach, and a numerical model fed with the output needs
straight inlet and outlet sections anyway. pyRiverBed therefore extends the
centerline at both ends with straight segments tangent to the first and last
existing segment. The extension length is $\lambda/10$ in Mode 1 and $B$ in
Mode 2. The extended polyline is then smoothed and resampled again.

---

## 4. Curvature

*Module: `pyriverbed.planform`*

For a centerline given as discrete points, pyRiverBed's default estimator is
the **`arctan2` method**: the local tangent direction is computed by forward
and backward differences,

```math
\theta^{+}_i = \mathrm{atan2}(y_{i+1}-y_i, x_{i+1}-x_i), \qquad
\theta^{-}_i = \mathrm{atan2}(y_i-y_{i-1}, x_i-x_{i-1})
```

and the curvature is the direction change per unit arc length, centred:

```math
C_i = \frac{2\left(\theta^{+}_i - \theta^{-}_i\right)}{s_{i+1} - s_{i-1}}
```

Two alternative estimators are implemented and available for comparison:

* **law of cosines** — the same idea, with the angles measured against a fixed
  reference direction;
* **three-point circumcircle** — the curvature is $1/R$ where $R$ is the
  radius of the circle through $P_{i-1}, P_i, P_{i+1}$,
  $R = abc/(4\sqrt{p(p-a)(p-b)(p-c)})$ with $a,b,c$ the triangle sides and
  $p$ the semi-perimeter. This one is unsigned, so it cannot distinguish a
  left turn from a right turn, which is why it is not the default.

**De-spiking.** After the estimate, any node whose curvature departs from the
mean of its neighbours by more than five times the neighbour-to-neighbour
difference is replaced by that mean:

```math
\text{if } \left|C_i - \frac{C_{i-1}+C_{i+1}}{2}\right| >
5\left|C_{i-1}-C_{i+1}\right| \Rightarrow 
C_i \leftarrow \frac{C_{i-1}+C_{i+1}}{2}
```

**Filtering.** The curvature signal is then filtered by a 5-point
second-order Savitzky–Golay filter followed by a 5-point first-order one (a
5-point moving average). Curvatures below a threshold $\varepsilon = 10^{-8}$
are snapped to zero so that straight reaches are exactly straight.

---

## 5. Curvature phase lag

*Module: `pyriverbed.planform`*

If the bed responded to the *local* curvature only, the deepest point of every
bend would sit exactly at the apex. It does not: in real rivers and in flume
experiments the pool is displaced downstream of the apex, typically by one to
several channel widths. The reason is that the secondary flow cell carries
*memory* — it takes a finite streamwise distance to develop after the flow
enters a bend, and it persists past the exit.

Linear bend theory makes this precise. The near-bank velocity perturbation
responds to the upstream curvature history through a convolution with an
exponentially decaying kernel,

```math
\mathcal{C}_{\text{lag}}(s) \propto 
\int_0^\infty C(s-\xi) e^{-a_2 \xi} \mathrm{d}\xi
```

pyRiverBed approximates this convolution with a discrete,
**linearly decaying, upstream-only weighted average** of the curvature:

```math
C^{\text{lag}}_i = \sum_{j=0}^{M-1} w_j C_{i-j},
\qquad
w_j = \frac{2}{M} - \frac{2j}{M(M-1)}
```

The weights are a normalised triangular kernel: they are largest at the local
node ($w_0 = 2/M$), decay linearly upstream, reach zero at $j = M-1$, and sum
to exactly one,

```math
\sum_{j=0}^{M-1} w_j = 2 - \frac{2}{M(M-1)}\cdot\frac{M(M-1)}{2} = 1
```

so the operator preserves the mean of the curvature signal and only shifts and
smooths it. The window length is set in channel widths by the **lag
strength** $\ell$:

```math
M = \left\lfloor \frac{\ell B}{\Delta s} \right\rfloor
```

so $\ell$ is directly interpretable: $\ell = 4$ means "the bed at this
location remembers the last four channel widths of curvature". Larger $\ell$
pushes the pools further downstream and damps their amplitude. Setting the lag
off recovers the purely local, apex-centred bed.

A triangular kernel rather than the theoretical exponential is a deliberate
simplification: it is compactly supported (no tail to truncate), exactly
normalised, and its first moment — which is what sets the phase shift — is
$\bar{j} = (M-1)/3$, i.e. one third of the window, a clean rule of thumb.

---

## 6. Bed topography

*Module: `pyriverbed.bed`*

### 6.1 Transverse bed slope

The transverse bed slope is taken proportional to the local curvature:

```math
S_T = A H C \xi_{S_T}
```

This is the classical closure of Ikeda et al. (1981) and Odgaard (1986):
the balance between the gravitational pull on a grain on a laterally sloping
bed and the transverse drag exerted on it by the secondary flow gives a
transverse slope proportional to $H/R = HC$. The constant of proportionality
$A$ is the **scour factor**, for which pyRiverBed uses Beck's (1988)
width-dependent form:

```math
A = 3.8\left[ 1 + \frac{\beta}{6.96} 
\exp\left(-\frac{6.96}{\beta}\right)\right],
\qquad \beta = \frac{b}{H}
```

$A \to 3.8$ for narrow, deep channels ($\beta \to 0$) and increases with
$\beta$: wide shallow channels scour their outer banks more aggressively
relative to their mean depth. For a laboratory flume with $\beta \approx 2$,
$A \approx 3.83$; for a large sand-bed river with $\beta \approx 10$,
$A \approx 6.5$.

$\xi_{S_T}$ is a purely empirical **transverse slope corrector** in $(0, 1]$.
It exists because the linear closure above systematically over-predicts the
transverse slope in sharp bends, where the response saturates. Reducing
$\xi_{S_T}$ scales the whole bar–pool relief down without touching its shape,
which is the knob to turn when calibrating against measured bathymetry.

**Sign convention.** With $n$ positive to the left of the flow and $C$ positive
for a left turn, this sign puts the pool against the **outer** bank, which is
where it belongs: a left-turning reach has its outer bank on the right, and
$S_T > 0$ makes the flow deepen towards $n < 0$.

> **Changed in v2.0.** v1.x carried the opposite sign here, which put the pool
> on the *inner* bank, and cancelled the error by shipping with its transverse
> flip switched on. v2 fixes the sign and defaults the flip to off, so the
> out-of-the-box result is the same one v1.x users actually got. A v1.x
> steering file is read with its `FLIPTRANS` value inverted, which reproduces
> v1.x output exactly. `flip.transverse` is now what its name suggests: a
> cosmetic mirror.

### 6.2 The Beck transverse profile

Given $S_T$, the **flow depth** across the section is

```math
\zeta(n) =
\begin{cases}
h_c - S_T n, & S_T n < 0 \quad \text{(pool side, outer bank)}\\[6pt]
h_c \exp\left(-\frac{S_T n}{H}\right), & S_T n > 0
\quad \text{(bar side, inner bank)}\\[6pt]
h_c, & n = 0
\end{cases}
```

This asymmetry is the heart of the Beck (1988) formulation and it is
physically motivated:

* Towards the **outer bank** the bed is a scour surface cut by a strong,
  quasi-uniformly loaded secondary current, and it comes out essentially
  **planar** — a straight transverse slope, which is exactly what surveyed
  pool cross-sections show.
* Towards the **inner bank** the bed is a *depositional* surface. The point
  bar is built by sediment carried up the slope, and as the water shallows the
  transport capacity falls off multiplicatively, not linearly. The result is
  the concave-up, **exponentially** shoaling bar face that is the signature
  shape of a point bar.

The two branches meet continuously at $n = 0$ with the common value $h_c$.

### 6.3 Centerline depth from mass conservation

$h_c$ is not a free parameter — it is fixed by requiring that the
cross-sectional area be preserved. The channel width is constant and the
reach-averaged depth is $H$, so every section must satisfy

```math
\int_{-b}^{b} \zeta(n) \mathrm{d}n = 2 b H
```

Substituting the two branches and solving for $h_c$ gives the closed form used
by the code:

```math
h_c = 
\frac{4 b H \vert S_T\vert - S_T^{2} b^{2}}
{2 b \vert S_T\vert + 2H - 
2H\exp\left(-\frac{\vert S_T\vert b}{H}\right)}
```

This is what guarantees that the synthetic bed is *self-consistent*: deepening
the pool automatically raises the bar by the compensating amount, so the
reach-averaged depth stays at $H$ and the generated bed can be handed to a
flow solver without re-balancing the discharge. Using $\vert S_T \vert$
makes the expression independent of the turning direction; the profile simply
mirrors.

In the straight-channel limit $S_T \to 0$ the expression is regular and gives
$h_c \to H$, i.e. a flat bed of uniform depth $H$, as it must. Numerically,
$\vert S_T \vert$ is floored at $\varepsilon = 10^{-8}$ to keep the
divisions safe.

### 6.4 Bed elevation

Flow depth becomes bed elevation by subtracting from the reference water
surface and adding back the longitudinal fall over the reach:

```math
z(s,n) = H - \zeta(s,n) + I\left(s_{\max} - s\right)
```

so $z = 0$ is the bed of an equivalent straight channel at the downstream end,
the bed rises upstream at the channel slope $I$, and pools are negative while
bars are positive. Setting $I = 0$ produces a horizontal-datum bed, which is
usually what you want for a flume design or for a purely geometric comparison.

Optionally the whole section can be **flipped transversely**, which mirrors the
bed about the centerline. This is the switch to use when your imported
centerline is digitised in the opposite direction from the one you intended,
or when you want the mirror image of a flume.

---

## 7. From 1D centerline to 2D channel

*Module: `pyriverbed.geometry`*

The bed is computed on a curvilinear $(s, n)$ grid, but the deliverable is a
point cloud in real $(x, y)$ coordinates. pyRiverBed builds it by **polyline
offsetting**: for each offset distance $L_k = k \Delta n$, $k = 1 \ldots N$,
it constructs the left and right polylines parallel to the centerline at
distance $L_k$, and attaches to their nodes the bed elevations from the
corresponding columns of the $(s,n)$ bed array.

Naive offsetting — displacing each node along its own normal — fails at the
sharp corners that a discretised centerline always has: the offset polyline
self-intersects on the inside of a bend and gaps open on the outside. The
implementation instead uses the **segment-intersection** method: each *segment*
is offset as a whole, and the offset node is placed at the intersection of the
two adjacent offset segments, found by solving the $2\times2$ linear system

```math
\begin{bmatrix}
\Delta y_{i-1} & -\Delta x_{i-1}\\
\Delta y_{i} & -\Delta x_{i}
\end{bmatrix}
\begin{bmatrix} x^{\text{off}}_i \\ y^{\text{off}}_i \end{bmatrix}
=
\begin{bmatrix}
\Delta y_{i-1} x^{0}_{i-1} - \Delta x_{i-1} y^{0}_{i-1}\\
\Delta y_{i} x^{0}_{i} - \Delta x_{i} y^{0}_{i}
\end{bmatrix}
```

where $(x^0, y^0)$ are the naively offset segment endpoints. Nearly collinear
segments (normalised dot product $> 1 - 10^{-10}$) would make the system
singular and are handled by falling back to the naive offset. This keeps the
channel width *exactly* constant, which is the defining assumption of the
whole framework.

The point cloud is assembled centerline first, then in pairs of left/right
offsets working outwards, giving the node ordering that the mesh generator
in §8 relies on. The outermost pair of offset polylines is the pair of
**banklines**, and it can be written out as a closed polygon.

---

## 8. Finite element mesh and boundary conditions

*Module: `pyriverbed.mesh`*

Because the point cloud is a structured
$(n_{\text{row}} \times n_{\text{col}})$ grid in disguise —
$n_{\text{row}}$ streamwise nodes,
$n_{\text{col}} = 2N+1$ transverse nodes — a triangulation can be written down
directly, with no Delaunay step. Each structured quadrilateral is split into
two triangles, giving

```math
n_{\text{ele}} = 2 (n_{\text{row}} - 1)(n_{\text{col}} - 1)
```

triangles whose connectivity follows from the node ordering of §7. Because the
streamwise and transverse spacings were made equal in §3.3, these triangles
are close to equilateral.

The following files are produced:

| File | Format | Use |
|:---|:---|:---|
| `*_topo.xyz` | 3-column ASCII point cloud | load in Blue Kenue, GIS, any interpolator |
| `*_boundary.i2s` | closed polyline | the banklines, for Blue Kenue |
| `*_mesh.t3s` | Blue Kenue T3 mesh | build a TELEMAC Selafin geometry file |
| `*_mesh.dat` | Tecplot FETRIANGLE | visualisation |
| `*_BC.cli` | TELEMAC boundary conditions | use directly in a TELEMAC run |
| `*_BC.bc2` | Blue Kenue BC | inspect and edit BC codes and metadata |

The boundary is walked in the order inlet, outlet, left bank, right bank. The
inlet section is tagged as a prescribed-discharge, free-surface boundary
(TELEMAC codes `4 5 5`) and the outlet as a prescribed-elevation boundary
(codes `5 4 4`); the banks stay as closed walls (codes `2 2 2`). Those are the
conventional choices for a steady subcritical run and are meant to be a
starting point, not a prescription.

---

## 9. Meander migration

*Module: `pyriverbed.migration`*

### 9.1 Linear bend theory

Meander growth is driven by the **excess velocity near the outer bank**.
Ikeda, Parker & Sawai (1981) linearised the shallow-water equations about a
straight-channel base state, closed the transverse slope as in §6.1, and
obtained a first-order ordinary differential equation for the dimensionless
near-bank velocity perturbation $u_b$ along the channel. Its solution is the
sum of

* a **boundary-condition transient** that decays exponentially downstream from
  the inlet,
* a term proportional to the **local curvature**, and
* a term proportional to the **convolution of the upstream curvature** with an
  exponential memory kernel — the same integral that §5 approximates.

pyRiverBed uses this closed form directly. Working in variables normalised by
the channel width ($\tilde{s} = s/B$, $\tilde{C} = CB$, and so on):

```math
u_b(\tilde{s}) = a_1 e^{-a_2 \tilde{s}} + a_3 \tilde{C}
 + a_4 \tilde{C}^{\text{lag}}
```

with

```math
\begin{aligned}
a_1 &= U_{b0} (2\xi - 1) + \chi C_0, \qquad \xi \sim \mathcal{U}(0,1)\\
a_2 &= 2 C_{f0} \beta \chi\\
a_3 &= -\chi\\
a_4 &= C_{f0} \beta\left[\chi^{5} Fr_0^{2}
 + (A+1) \chi^{2}
 + 5\sqrt{C_{f0}}\left(A + \chi^{2} Fr_0^{2}\right)\right]
\end{aligned}
```

where

```math
\chi = \left(\frac{\Vert P_L - P_0\Vert}{L}\right)^{1/3}
 = \Omega^{-1/3}
```

$\chi$ is the feedback of the planform on the hydraulics. As the river
meanders its length $L$ grows while its valley length stays fixed, so its
slope falls by the factor $\Omega$; through a Chézy-type normal-flow relation
the velocity scale falls as the cube root of the slope, hence
$\chi = \Omega^{-1/3}$. This is what makes long simulations self-limiting:
a very sinuous river migrates more slowly, all else equal.

**The sign structure matters.** The local-curvature coefficient $a_3$ is
negative while the lagged coefficient $a_4$ is positive and, for realistic
parameters, several times larger — with $C_{f0} = 0.01$, $\beta = 21$,
$Fr_0 = 0.1$ one gets $a_3 = -1$ against $a_4 \approx 4$. It is the
*phase-shifted* term that dominates, and that is precisely why meanders both
**grow** in amplitude and **translate downstream** rather than simply
oscillating in place. A model without the lag would not meander.

$U_{b0}$ and $C_0$ set the inlet perturbation: $C_0$ a deterministic one and
$U_{b0}$ the amplitude of a random one, redrawn every time step. A
straight channel is a fixed point of the deterministic model, so
$U_{b0} > 0$ is what lets meanders grow out of an initially straight
channel (`straight.txt`) — the numerical analogue of the natural noise that
seeds a bend instability.

### 9.2 Bank migration

Bank retreat is taken proportional to the near-bank excess velocity — the
classical Ikeda et al. (1981) erosion law, whose coefficient was calibrated
across many rivers by Hasegawa and by Pizzuto. Every node is displaced along
the local normal:

```math
\tilde{x} \leftarrow \tilde{x} + E_0 u_b \Delta t \sin\theta,
\qquad
\tilde{y} \leftarrow \tilde{y} - E_0 u_b \Delta t \cos\theta
```

$(\sin\theta, -\cos\theta)$ is the unit normal to the right of the flow
direction, so a positive $u_b$ moves the centerline to the right.

Note the model is *kinematic in width*: both banks move together and the
channel width never changes. Width adjustment, floodplain heterogeneity and
bank stratigraphy are outside the scope of the framework.

After displacement the centerline is re-smoothed and resampled (§3.3, §3.4),
which both removes the small-scale roughness that the pointwise displacement
introduces and keeps the node spacing uniform as the reach lengthens. The
number of smoothing passes per step matters more than it looks: too few and
node-scale curvature spikes survive and inflate the sinuosity with wiggles that
are not bends; too many and the smoother diffuses the meanders faster than the
model grows them, so the reach decays towards a straight line. Two passes is
the knee, and results are insensitive up to roughly eight.

### 9.3 The upstream boundary

The linearised theory is an *initial value problem in $s$*: $u_b$ at any point
depends on the curvature upstream of it. At the first node there is no upstream
reach, so the inlet perturbation $a_1$ enters as a boundary value and the
solution near the inlet is whatever that boundary value dictates.

Left free, this is unstable. The perturbation is nearly uniform over the first
$1/a_2 \approx 2$–3 channel widths, and any small bend it produces there is
amplified by the $a_4$ term without an upstream reach to spread it over. The
inlet develops a hook whose curvature quickly exceeds anything the linear
theory can describe, and the hook then folds onto the reach and triggers
spurious cutoffs.

pyRiverBed therefore fixes the end nodes, which is the standard treatment of a
free boundary in a meander migration model. The displacement is multiplied by a
raised-cosine taper

```math
w(s) = \frac{1}{2}\left[1 - \cos\left(\pi \min\left(
\frac{\min(s, L-s)}{\Lambda_{\text{taper}} B}, 1\right)\right)\right]
```

which is zero at both ends and reaches one over `end_taper_widths` channel
widths. A raised cosine rather than a straight ramp so that $w$ and $w'$ are
both continuous where the taper meets the interior; a linear ramp would leave a
weak kink there. Setting `end_taper_widths = 0` recovers the free ends, which
are only stable under heavy per-step smoothing.

The interior physics is untouched: the taper is exactly one wherever
$\min(s, L-s) \geq \Lambda_{\text{taper}} B$. The cost is that the reach cannot
translate at its very ends, so the first and last couple of widths should be
treated as buffer and not interpreted.

### 9.4 Diagnostics

Two time series are written for the whole run:

* **sinuosity** $\Omega_t = L_t / \Vert P_L - P_0 \Vert_t$ — the standard
  measure of meander development, which grows, then plateaus as cutoffs start
  to shorten the channel as fast as migration lengthens it;
* **mean migration rate** — the mean nodal displacement per time step,
  $\frac{1}{n}\sum_i \Vert P_i^{t+1} - P_i^{t}\Vert$, useful for checking
  that $E_0 \Delta t$ gives a physically sensible rate and that the run is not
  taking steps so large that the planform is unstable.

---

## 10. Cutoffs

*Module: `pyriverbed.migration`*

A meander cannot grow forever. Sooner or later the channel short-circuits, the
bypassed loop is abandoned and becomes an **oxbow lake**, and the channel
length — and hence the sinuosity — drops abruptly. pyRiverBed models both
mechanisms observed in the field, and they are genuinely different processes.

### 10.1 Neck cutoff

A neck cutoff happens when a meander loop grows until its two limbs *touch*.
It is a purely geometric event, so pyRiverBed detects it geometrically:
a cutoff is declared at the first pair of nodes $(i, j)$ that are

* **far apart along the channel**: $j - i > 4N$ nodes, i.e. more than
  $2B$ of arc length, and
* **close together in space**: $\Vert P_i - P_j \Vert < B$, i.e. less than
  one channel width.

The reach between them is removed and returned as an oxbow lake, and the arc
length is recomputed. No stochastic ingredient is involved — a neck cutoff is
an inevitable consequence of sustained growth.

`end_margin_widths` of each end is excluded from the search, for the same
reason it is excluded for chute cutoffs (§10.2): those stretches carry the
artificial straight extensions of §3.5 and the inlet boundary treatment of
§9.3, so a qualifying pair found there reflects the boundary condition rather
than a meander loop closing on itself.

### 10.2 Chute cutoff

A chute cutoff is different in kind. The channel does not close on itself;
instead the flow carves a **new, shorter channel across the floodplain or
across the point bar**, abandoning the bend while its limbs are still well
separated. Reviews of the process (Constantine et al., 2010; Grenfell et al.,
2012; van Dijk et al., 2012; Zinger et al., 2011) identify several routes to
it — headward incision of a swale on the bar surface, headcutting back from
the downstream limb, or a mid-channel bar diverting flow into an inner
branch — and they share two features that matter for a model:

1. **A chute needs a slope advantage.** The new path is taken only if it is
   appreciably steeper than the along-channel path, which means it must be
   reasonably well aligned with the *down-valley* direction and must bypass a
   substantial length of channel.
2. **A chute needs a trigger.** Unlike a neck cutoff, chute initiation is
   contingent on events the model does not resolve — a flood of the right
   magnitude, a weak spot in the bank, a gap in the riparian vegetation. It is
   properly treated as a **stochastic** process.

pyRiverBed therefore models chute cutoffs as *conditionally random*: the
geometry decides where a chute is **possible**, and a random draw decides
whether one **happens**.

**Where a chute is possible.** The filtered curvature signal gives the bend
structure of the reach. Sign changes of $C$ are the **inflection points** (the
crossings); the point of locally maximum $\vert C \vert$ between two
consecutive inflection points is a **bend apex**. Either family can serve as
the entrance and exit of a chute, selected by `entrance`:

* `apex` — the chute leaves and rejoins the channel at bend apexes. This is
  the bar-surface / swale-incision route.
* `inflection` — the chute leaves and rejoins at crossings, cutting a whole
  bend out of the planform.

A candidate chute connects entrance point $k$ to entrance point $k +$ `span`,
so `span = 2` bypasses one full meander loop, `span = 4` bypasses two, and so
on.

**The geometric criteria.** A candidate is admissible if all four hold:

```math
\underbrace{j - i > \Lambda_{\min}\cdot 2N}_{\text{long enough}}
\qquad
\underbrace{\frac{s_j - s_i}{\Vert P_j - P_i\Vert} \geq \Omega_{\min}}
_{\text{worth taking}}
\qquad
\underbrace{\alpha(i,j) \leq \alpha_{\max}}_{\text{aligned enough}}
\qquad
\underbrace{\mu \leq k \leq n_{\text{ent}} - \text{span} - \mu}
_{\text{away from the ends}}
```

**Long enough.** $\Lambda_{\min}$ (`min_length_widths`) is the minimum length of
the bypassed reach in channel widths. Since the node spacing is
$\Delta n = b/N$, a reach of $j-i$ nodes is $(j-i) B/(2N)$ long, so the
node-count threshold is $\Lambda_{\min} \cdot 2N$. It encodes "a chute is only
worth taking if it saves a meaningful amount of distance".

**Worth taking.** $\Omega_{\min}$ (`min_sinuosity`) is the minimum sinuosity of
the bypassed reach — the ratio of the channel distance from $i$ to $j$ to the
straight chute chord across it. **This is the slope advantage itself**, and
therefore the most direct statement of why a chute forms at all: the new path is
taken only because it is appreciably steeper than the old one. Without this test
a chute can be carved across a reach that is already straight, gaining the flow
nothing, which is not something rivers do.

**Aligned enough.** $\alpha(i,j)$ (`max_valley_angle`) is the angle between the
**chute chord** and the **valley axis**. Writing $\mathbf{c} = P_j - P_i$ for
the chord and $\mathbf{v} = P_L - P_0$ for the valley axis,

```math
\alpha = \arctan\frac{\vert \mathbf{c}\times\mathbf{v}\vert}
{\vert \mathbf{c}\cdot\mathbf{v}\vert}
 \in [0^\circ, 90^\circ]
```

A small $\alpha$ means the chute runs straight down the valley and therefore
enjoys the largest possible slope advantage; a chute perpendicular to the valley
has none and is rejected. Because the valley axis is measured from the
centerline itself, this test is invariant to how the reach is rotated in the
coordinate system.

**Away from the ends.** $\mu$ (`end_margin`) keeps chutes away from the inlet
and outlet, where the planform is contaminated by the straight extensions of
§3.5 and by the inlet transient of §9.1.

**Whether a chute happens.** At every time step after a spin-up of
`start_step` steps, a chute cutoff is triggered with probability
`frequency`. `frequency = 0.1` means a 10% chance per time step, i.e. a mean
recurrence of ten time steps. When triggered, **all** admissible candidates are
collected and one is drawn uniformly at random and carved; if there are none,
nothing happens.

Note that the realised cutoff rate is bounded from above by the *supply* of
admissible bends, not only by `frequency`: every cutoff straightens the reach,
and once no bend still satisfies the criteria, raising `frequency` further
changes nothing until migration has rebuilt the meander belt. Section 9 of the
demo notebook shows this saturation directly.

> **Changed in v1.1.0 and v2.0.** v1.x drew a *single* candidate and tested it,
> so a step in which the drawn bend happened to fail produced no cutoff even
> when other bends would have qualified. It also computed the valley angle with
> a call that always raised inside a bare `except`, so before v1.1.0 the
> alignment criterion had no effect whatsoever, and every one of these
> thresholds was a literal constant in the source rather than an input.

The spin-up exists because a nearly straight initial channel has no meaningful
bend structure — allowing cutoffs before the bends have grown would just
truncate the reach at random.

**Carving.** Identical bookkeeping to a neck cutoff: nodes $i+1 \ldots j-1$
become the oxbow lake, the centerline becomes $P_0 \ldots P_i, P_j \ldots P_L$,
and the arc length is recomputed. The kink at the junction is removed by the
smoothing and resampling that follow (§3.3, §3.4) — which is also, physically,
the widening and reworking of the young chute channel.

At most one cutoff of each kind is carved per time step, and a chute cutoff is
skipped in any step where a neck cutoff has already fired, since the node
indices the chute search returned refer to the pre-cutoff centerline.

---

## 11. Assumptions and limits of validity

pyRiverBed is a *geometric–analytical* generator, not a morphodynamic solver.
Being explicit about what it does not do is the best guide to when it can be
trusted.

**Constant width.** The width is a fixed input everywhere and at all times.
Real meanders widen at bends and narrow at crossings, and chute cutoffs in
particular are strongly modulated by width variation. If width variation is
central to your question, this is the wrong tool.

**Local, linear curvature–bed coupling.** The bed at a section depends only on
the (phase-lagged) curvature at that section. There is no sediment continuity
equation, no bar migration, no grain-size sorting, and no upstream sediment
supply. The consequence is that the bed is an **equilibrium** bed: it is the
topography the reach would have if the current planform persisted long enough,
not a transient state.

**Linearised hydrodynamics.** The migration model is a first-order theory
about a straight base state. It is quantitatively reliable for mild curvature
($CB$ of order 0.3 or less) and moderate width-to-depth ratios. For very sharp bends
the linear closure over-predicts the transverse slope — which is exactly what
$\xi_{S_T}$ is there to absorb — and the near-bank velocity is no longer a
sufficient description of the flow field.

**No floodplain.** There is no floodplain topography, no stratigraphy, no
erodibility field, no vegetation. Oxbow lakes are recorded for plotting but do
not influence subsequent migration, and abandoned bends do not resist being
re-occupied.

**Stochastic chute cutoffs are a parameterisation, not a mechanism.** The
`frequency` parameter stands in for flood history, bank strength and bar
topography — everything the model does not resolve. It should be *calibrated*
against an observed cutoff rate for the river of interest, and results should
be interpreted as one realisation of an ensemble, not as a prediction. Run the
model several times to see the spread.

**Rigid time stepping.** There is no stability control on $E_0\Delta t$. If
the per-step displacement approaches the node spacing, the planform will
develop numerical wiggles that the smoother then has to fight. Watch the mean
migration rate diagnostic.

---

## 12. Choosing parameters

| Parameter | Typical range | How to choose it |
|:---|:---|:---|
| `width`, `depth` | – | Measured, or bankfull estimates. Their ratio $\beta = B/2H$ drives the scour factor. |
| `slope` | 0 – 10⁻³ | Set to 0 for a horizontal datum; otherwise the reach-averaged channel slope. |
| `transverse_slope_corrector` | 0.5 – 1.0 | Start at 1.0. Lower it if the synthetic pools are deeper than the surveyed ones — this is the main bar–pool relief calibration knob. |
| `ds` | $\lambda/500$ – $\lambda/200$ | Mode 1 only. Fine enough to resolve the bend apexes. |
| `n_offsets` | 5 – 20 | Sets the transverse resolution *and*, through §3.3, the streamwise one. 10 is a good default; raise it for a mesh you intend to run a solver on. |
| `smoothing_level` | 0 – 100 | Mode 2 only, and the single most important number for a real centerline. Raise it until the curvature signal looks like a bend sequence rather than noise, and no further. |
| `lag_strength` | 2 – 6 | In channel widths. Calibrate against the observed offset between bend apex and deepest point; 4 is a reasonable default. |
| `cf0` | 0.01 – 0.03 | Reach-averaged friction coefficient. |
| `fr0` | 0.1 – 0.5 | Reach-averaged Froude number; subcritical. |
| `e0` | 10⁻⁸ – 10⁻⁵ s⁻¹ | Bank erosion coefficient. Calibrate against an observed migration rate — this is the *only* rate-setting parameter, so tune it and leave the rest. |
| `dt` | 3600 – 86400 s | Time step. Only the product $E_0\Delta t$ matters for the outcome; $\Delta t$ alone sets the time labelling of the output. |
| `ub0` | 0 – 5 | Inlet noise amplitude. Must be non-zero to grow meanders from a straight channel. |
| `migration_smoothing_level` | 2 – 8 | Passes per time step. Below 2, node-scale curvature spikes survive; above ~8 the smoother diffuses the meanders. Not the same knob as `smoothing_level`. |
| `end_taper_widths` | 2 | Length of the fixed-end ramp, in channel widths, covering the inlet transient $1/a_2$. 0 gives free ends, stable only under heavy per-step smoothing. |
| `end_margin_widths` | 2 | Reach ends excluded from neck cutoff detection. Match it to `end_taper_widths`. |
| `chute.frequency` | 0.001 – 0.5 | Per-time-step probability. Calibrate against an observed cutoff recurrence: `frequency` $\approx \Delta t / T_{\text{recurrence}}$. |
| `chute.start_step` | 10 – 50 % of `n_steps` | Long enough for bends to develop. |
| `chute.max_angle` | 15° – 45° | Lower is more selective. Set to 90° to disable the alignment test. |
| `chute.min_length` | 5 – 20 | Bypassed reach length in channel widths. |
| `chute.span` | 2 | 2 bypasses one meander loop. Use 4 for the occasional double-loop cutoff. |

---

## 13. References

Abad, J. D., & Garcia, M. H. (2009). Experiments in a high-amplitude Kinoshita
meandering channel: 1. Implications of bend orientation on mean and turbulent
flow structure. *Water Resources Research*, 45(2), W02401.
doi:[10.1029/2008WR007016](https://doi.org/10.1029/2008WR007016)

Beck, S. (1988). Computer-simulated deformation of meandering river channels.
Ph.D. thesis / technical report.

Constantine, J. A., McLean, S. R., & Dunne, T. (2010). A mechanism of chute
cutoff along large meandering rivers with uniform floodplain topography.
*GSA Bulletin*, 122(5–6), 855–869.
doi:[10.1130/B26560.1](https://doi.org/10.1130/B26560.1)

Grenfell, M., Aalto, R., & Nicholas, A. (2012). Chute channel dynamics in
large, sand-bed meandering rivers. *Earth Surface Processes and Landforms*,
37(3), 315–331. doi:[10.1002/esp.2257](https://doi.org/10.1002/esp.2257)

Ikeda, S., Parker, G., & Sawai, K. (1981). Bend theory of river meanders.
Part 1. Linear development. *Journal of Fluid Mechanics*, 112, 363–377.
doi:[10.1017/S0022112081000451](https://doi.org/10.1017/S0022112081000451)

Kinoshita, R. (1961). *Investigation of channel deformation in Ishikari
River*. Report, Bureau of Resources, Department of Science and Technology,
Japan.

Langbein, W. B., & Leopold, L. B. (1966). River meanders — theory of minimum
variance. *USGS Professional Paper* 422-H.

Li, Z., & Garcia, M. H. (2021). pyRiverBed: A Python framework to generate
synthetic riverbed topography for constant-width meandering rivers.
*Computers & Geosciences*, 152, 104755.
doi:[10.1016/j.cageo.2021.104755](https://www.doi.org/10.1016/j.cageo.2021.104755)

Motta, D., Abad, J. D., Langendoen, E. J., & Garcia, M. H. (2012). A
simplified 2D model for meander migration with physically-based bank
evolution. *Geomorphology*, 163–164, 10–25.
doi:[10.1016/j.geomorph.2011.06.036](https://doi.org/10.1016/j.geomorph.2011.06.036)

Odgaard, A. J. (1986). Meander flow model. I: Development. *Journal of
Hydraulic Engineering*, 112(12), 1117–1136.
doi:[10.1061/(ASCE)0733-9429(1986)112:12(1117)](https://doi.org/10.1061/(ASCE)0733-9429(1986)112:12(1117))

Parker, G., Sawai, K., & Ikeda, S. (1982). Bend theory of river meanders.
Part 2. Nonlinear deformation of finite-amplitude bends. *Journal of Fluid
Mechanics*, 115, 303–314.
doi:[10.1017/S0022112082000767](https://doi.org/10.1017/S0022112082000767)

Rowley, T., Konsoer, K., Langendoen, E. J., Li, Z., Ursic, M., & Garcia, M. H.
(2021). Relationship of point bar morphology to channel curvature and planform
evolution. *Geomorphology*, 375, 107541.
doi:[10.1016/j.geomorph.2020.107541](https://www.doi.org/10.1016/j.geomorph.2020.107541)

Savitzky, A., & Golay, M. J. E. (1964). Smoothing and differentiation of data
by simplified least squares procedures. *Analytical Chemistry*, 36(8),
1627–1639. doi:[10.1021/ac60214a047](https://doi.org/10.1021/ac60214a047)

van Dijk, W. M., van de Lageweg, W. I., & Kleinhans, M. G. (2012).
Experimental meandering river with chute cutoffs. *Journal of Geophysical
Research: Earth Surface*, 117, F03023.
doi:[10.1029/2011JF002314](https://doi.org/10.1029/2011JF002314)

Zinger, J. A., Rhoads, B. L., & Best, J. L. (2011). Extreme sediment pulses
generated by bend cutoffs along a large meandering river. *Nature Geoscience*,
4, 675–678. doi:[10.1038/ngeo1260](https://doi.org/10.1038/ngeo1260)
