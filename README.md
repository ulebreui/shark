# SHARK

SHARK is a 2D finite-volume code for gas and multi-species dust dynamics. It is designed for
protoplanetary disks and local dust–gas problems. It solves the hydrodynamics of the gas and of N
pressureless dust fluids coupled by drag, with optional dust growth, in Cartesian or polar (face-on disk)
geometry. It is parallelised with OpenMP.

## Features

**Gas**
- Second-order MUSCL-Hancock scheme (unsplit predictor + Godunov corrector)
- Riemann solvers: LLF, HLL, HLLC
- Slope limiters: minmod, van Leer, MC
- Isothermal or adiabatic equation of state

**Dust**
- N dust species, treated as pressureless fluids
- Dust Riemann solvers: Huang & Bai (default), LLF, HLL, HLLC
- Drag with back-reaction, integrated implicitly (Krapp et al. scheme)
- Dust growth:
  - full Smoluchowski coagulation/fragmentation solver
  - single-size (Stepinski-like) growth carried by a passive scalar, capped at the fragmentation-limited size

**Disk geometry** (polar, face-on)
- Linear, logarithmic or lin–log radial grid
- FARGO orbital advection (Masset 2000): mean azimuthal velocity removed from the transport, followed by a conservative ring shift
- Planets on fixed circular orbits:
  - Plummer softening, mass ramp-up and indirect term
  - torques on each planet written to a file
- Wave-damping zones at the radial boundaries (de Val-Borro et al. 2006)

**Local setups**
- Shearing box for the streaming instability
- Vertical dust settling
- 1D validation problems

**Time step:** unsplit CFL condition, with a force constraint and a shear constraint when FARGO is used.

## Repository layout

```
bin/            Makefile (compile here)
src/
  shark.f90     main program and time loop
  commons/      parameters, allocation, namelist reading, outputs
  grid/         grid construction
  hydro/        hydro arrays
  solver/       predictor, Riemann solvers, Courant condition, source terms, FARGO
  dust/         dust arrays, drag, growth
  setups/       one folder per problem (setup.f90, setup_commons.f90, boundaries.f90, units.f90)
    disk_polar/ global disk (planets, FARGO, damping zones)
    SI/         streaming instability (shearing box)
    settling/   vertical dust settling
    tests/      1D validation problems
test/
  run_tests.py  analytical validation suite
  test-suite/   namelists of the test problems
```

## Compilation

SHARK needs a Fortran compiler (gfortran or ifort). All options are passed to `make` from `bin/`:

```bash
cd bin
make clean && make SETUP=disk_polar GEOM=2 GRIDSPACE=1 NX=256 NY=256 NDUST=1 NDUSTPSCAL=1 OPENMP=1 SOLVER=1 SOLVERDUST=0
```

| Option | Meaning |
|---|---|
| `SETUP` | problem folder in `src/setups/` (`disk_polar`, `SI`, `settling`, `tests`) |
| `GEOM` | `0` Cartesian, `2` polar disk |
| `GRIDSPACE` | radial grid for `GEOM=2`: `0` linear, `1` logarithmic, `2` lin–log |
| `NX`, `NY` | number of active cells (`NY=1` for 1D) |
| `NGHOST` | number of ghost cells (default 2) |
| `NDUST` | number of dust species (`0` = gas only) |
| `NDUSTPSCAL` | passive scalars per dust species (the first one is the grain size for Stepinski growth) |
| `SOLVER` | gas solver: `0` LLF, `1` HLL, `2` HLLC |
| `SOLVERDUST` | dust solver: `0` Huang & Bai, `1` LLF, `2` HLL, `3` HLLC |
| `OPENMP` | `1` to compile with OpenMP |
| `DEBUG` | `1` for bounds checking and floating-point traps |
| `COMP` | `"GNU"` (default) or `"INTEL"` |

Always run `make clean` when changing `SETUP` or any compile-time option.

## Running

```bash
export OMP_NUM_THREADS=16
cd my_run/
/path/to/bin/shark my_run.nml
```

The namelist file contains the groups `&GRID_PARAMS`, `&SETUP_PARAMS` and, when `NDUST>0`, `&DUST_PARAMS`.
Examples are in `test/test-suite/`.

**Main `&GRID_PARAMS` entries**

| Entry | Meaning |
|---|---|
| `tend` | final time (in years for `disk_polar`, in code units otherwise) |
| `CFL` | Courant number |
| `freq_out` | output frequency, in number of steps |
| `slope_type` | limiter: `1` minmod, `2` van Leer, `3` MC |
| `gamma` | adiabatic index |
| `force_kick` | apply the external forces |
| `fargo` | FARGO orbital advection (`GEOM=2`) |
| `nrestart` | restart from this output |

**Main `&SETUP_PARAMS` entries for `disk_polar`**

| Entry | Meaning |
|---|---|
| `box_l` | outer radius (au) |
| `smooth_r` | inner radius (au) |
| `sigma_R0`, `R0_disk`, `R_out_disk` | surface density profile: `Σ = Σ0 (R0/r) exp(-r/R_out)` |
| `hoverr` | aspect ratio |
| `Mstar` | stellar mass (Msun) |
| `damping`, `damp_in_frac`, `damp_out_frac`, `damp_tau` | wave-damping zones |
| `nplanet`, `Mplanet(:)`, `Rplanet(:)`, `phi_planet(:)` | planets: number, masses (Msun), radii (au), initial azimuths |
| `nramp_planet(:)` | mass ramp-up duration, in planet orbits |
| `a_smooth_pl` | softening, as a fraction of the local scale height (0.4–0.6 is typical) |
| `indirect_term` | include the indirect term |
| `freq_torque` | number of steps between two torque outputs (`planet_torques.dat`) |

## Outputs

Each output is written to `output_XXXXX/`, as raw binary arrays (float64, stream access) over the active
cells, with `iy` as the outer loop:

| File | Contents |
|---|---|
| `info.dat` | grid size, time, units |
| `x`, `y` | positions (`x` is the radius in polar geometry) |
| `rho`, `v`, `vy`, `vz`, `P` | gas density (surface density in polar geometry), velocities, pressure |
| `rhod`, `vd`, `vdy`, `vdz` | dust density and velocities, all species one after the other |
| `dustpscal` | passive scalars (grain size for Stepinski growth) |
| `metrics.dat` | time spent in each part of the step |
| `planet_torques.dat` | per line: time, planet index, current mass, gas torque, dust torque (code units) |

Reading an output with NumPy:

```python
import numpy as np
NX, NY = 256, 256
rho = np.fromfile("output_00010/rho", dtype=np.float64).reshape(NY, NX)   # [iy, ix]
```

## Validation

`test/run_tests.py` builds SHARK, runs a set of problems with known solutions, compares them and prints
PASS/FAIL. It saves one figure per test in `test/validation_runs/figures/`.

```bash
cd test
python3 run_tests.py                 # standard tests
python3 run_tests.py --long          # + planet type-I torque
python3 run_tests.py --only sod dustywave
python3 run_tests.py --threads 8
```

| Test | Reference |
|---|---|
| `sod` | Sod shock tube vs exact Riemann solution |
| `soundwave` | isothermal sound wave, second-order convergence |
| `dustybox` | drag relaxation, 1 and 3 species (Laibe & Price 2011) |
| `dustywave` | linear gas + dust wave (Laibe & Price 2011) |
| `dustyshock` | stationary isothermal dusty shock (Lorén-Aguilar & Bate 2014) |
| `disk_ic_drift` | initial dust radial drift in a disk |
| `disk_equilibrium` | disk equilibrium, second-order convergence |
| `dust_drift_nsh` | dust–gas radial drift (Nakagawa et al. 1986) |
| `fargo_bump` | dust bump advected by the rotation, with and without FARGO |
| `si_linA` | streaming instability growth rate (Youdin & Johansen 2007, linA) |
| `planet_torque` | type-I torque (Paardekooper et al. 2010), with `--long` |

The suite requires Python 3 with NumPy. Matplotlib is optional and is only used for the figures.
The exit code is the number of failed tests, so the suite can be run in continuous integration.

## Adding a setup

Copy an existing folder of `src/setups/` and adapt these four files:

| File | Role |
|---|---|
| `setup_commons.f90` | module `setup_parameters` (setup parameters) |
| `setup.f90` | initial conditions, plus the routines called by the core: `read_setup_params`, `flag_continue`, `check_output`, `setup_preloop`, `setup_inloop`, `update_force_setup`, `compute_tstop` (if `NDUST>0`), `write_setup_info`, `read_setup_info` |
| `boundaries.f90` | `apply_boundaries` |
| `units.f90` | module `units` |

Then compile with `SETUP=<folder>`.

Fortran ignores case, so a local variable named like a module variable (`ip` and `iP`, `p_l` and `P_l`)
causes a compile error or, worse, silently overwrites a global index. Prefix setup-specific names
(e.g. `sod_rho_l`).
