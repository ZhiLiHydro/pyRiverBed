# Example input files

Run any of these from the repository root:

```bash
pyriverbed run examples/01_kinoshita_flume.ini -o out
```

| File | What it does |
|:---|:---|
| `01_kinoshita_flume.ini` | The Abad & Garcia (2009) laboratory flume: a 0.6 m wide Kinoshita channel. These are pyRiverBed's defaults, and it is the fastest way to check an installation. |
| `02_jurua_river.ini` | The Juruá River reach in Brazil, read from `example_centerlines/jurua.txt`. Shows Mode 2 and the effect of `smoothing_level` on a hand-digitised centerline. |
| `03_meander_migration_cutoffs.ini` | A 20 000-step migration run from a Kinoshita curve, with both neck and chute cutoffs. Produces the meander belt figures and animations. |
| `04_straight_to_meandering.ini` | Grows meanders out of an initially straight 210 m wide channel, driven by the random inlet perturbation `ub0`. The slowest of the five, and the source of the meander belt and art prints in the main README. |
| `05_chute_frequency.ini` | The infrequent-chute case of the three-way comparison in the main README. Set `chute_cutoff.enabled = no` for the neck-only case and `frequency = 0.01` for the frequent-chute case; the seed is fixed so the three are directly comparable. |

A few things worth knowing before you edit them:

* **`migration.e0 * migration.dt` is the displacement per time step, in channel
  widths.** It is the only rate-setting quantity. Keep it well below 0.1.
* **`curvature.migration_smoothing_level` belongs in 2–8.** Smoothing diffuses
  meanders, so a large value reapplied every step destroys bends faster than the
  migration model grows them; a single pass, on the other hand, leaves
  node-scale curvature spikes that inflate the sinuosity with wiggles that are
  not bends.
* **`chute_cutoff.frequency` is a probability per time step**, so calibrate it
  as `dt / recurrence_interval` for your river. With `dt` of one day, the
  `0.001` these files use is a chute cutoff every ~2.7 years.
* **`migration.end_taper_widths` and `neck_cutoff.end_margin_widths` keep the
  reach ends out of the physics.** Treat the first and last couple of channel
  widths as buffer and do not interpret them.
* **Set `migration.seed`** to make a stochastic run reproducible, and run several
  seeds before drawing conclusions from one.

See [`../THEORY_GUIDE.md`](../THEORY_GUIDE.md) for what every parameter means
physically, and `pyriverbed init` for a fully commented template.
