# Benchmark: 3D Advection-Reaction

This benchmark problem implements a 3D advection-reaction equation using the
RAJA performance portability layer with serial, CUDA, or HIP backends.

## Problem description

This code simulates the advection and reaction of three chemical species where
the reaction mechanism is a variation of the Brusselator problem from chemical
kinetics. The PDE system is given by
```
    u_t = -c grad(u) + A - (w+1) * u + v * u^2
    v_t = -c grad(v) + w * u - v * u^2
    w_t = -c grad(w) + (B - w) / epsilon - w * u
```
where `u`, `v`, and `w` are chemical concentrations, `c` is the advection speed,
`A` and `B` are the concentrations of chemical species that remain constant over
space and time, and `epsilon` is a parameter that varies the stiffness of the
system. The problem is solved on the domain `(x,y,z) = X` in `[0, X_max]^3`,
for times `t` in `[0,t_f]`. The initial condition is
```
    u(0,X) = A + p(X)
    v(0,X) = B / A + p(X)
    w(0,X) = 3.0 + p(X)
```
where the perturbation function is
```
    p(X) = alpha * e^( -((X-mu)^T sigma^{-1} (X-mu)) / (2*sqrt(|sigma| 8 pi^3)) )
```
with `alpha = 0.1`, `mu = 0.5 X_max`, and `sigma` is a diagonal matrix with
entries `0.25 X_max`.

Spatial derivatives are discretized with first-order upwind finite differences
on a uniform spatial grid. The system can be evolved in time using explicit,
implicit, or IMEX methods from ARKODE, Adams or BDF methods from CVODE, or BDF
methods from IDA. When using an IMEX method, advection is treated explicitly and
reactions implicitly.

The nonlinear system(s) that arise in each time step may be solved using a
global Newton method with a matrix-free GMRES linear solver or an Anderson
accelerated fixed-point method. When using an IMEX method, a custom task-local
nonlinear solver that leverages the locality of the reaction systems may also be
used.

## Options

Several command line options are available to change the problem parameters
as well as the integrator and solver options. A summary of the options are
listed below.

| Option                      | Description                                                                   | Default     |
|:----------------------------|:------------------------------------------------------------------------------|:------------|
| `--help`                    | Print the command line options and description                                | --          |
| `--dont-save`               | Do not save the solution to the disk                                          | Save        |
| `--output-dir <dir>`        | Directory where all output files will be written                              | `.`         |
| `--nout <int>`              | Number of output times                                                        | 40          |
| `--npts <int>`              | Number of mesh points in each direction                                       | 100         |
| `--npxyz <int> <int> <int>` | Number of MPI tasks in each direction (0 forces MPI to decide)                | 0 0 0       |
| `--xmax <realtype>`         | Maximum value of `x`, `y`, and `z` in :math:`X_max`                           | 1.0         |
| `--A <realtype>`            | Constant concentration of species `A`                                         | 1.0         |
| `--B <realtype>`            | Constant concentration of species `B`                                         | 3.5         |
| `--c <realtype>`            | Advection speed `c`                                                           | 0.01        |
| `--order <int>`             | Integration method order                                                      | 3           |
| `--method <method>`         | Integrator to use: `ERK`, `ARK-DIRK`, `ARK-IMEX`, `CV-BDF`, `CV-ADAMS`, `IDA` | `ARK-DIRK`  |
| `--nls <method>`            | Nonlinear Solver Method: `newton`, `tl-newton`, `fixedpoint`, `none`          | `newton`    |
| `--fpaccel <int>`           | Number of fixed point acceleration vectors                                    | 3           |
| `--nopre`                   | Disable preconditioning                                                       | False       | 
| `--fused`                   | Enabled fused operations                                                      | Off         |
| `--tf <realtype>`           | Final integration time `t_f`                                                  | 10.0        |
| `--rtol <realtype>`         | Relative tolerance                                                            | 1.0e-6      |
| `--atol <realtype>`         | Absolute tolerance                                                            | 1.0e-9      |

## Building and Running

To build the benchmark executables SUNDIALS must be configured with ARKODE,
CVODE, and IDA enabled and with MPI and RAJA support on. Additionally, either
CUDA or HIP support must be on to build executables utilizing NVIDIA or AMD
GPUs. See the installation guide for more details on configuring, building,
and installing SUNDIALS.

Based on the configuration the following executables will be built and installed
in the `<install prefix>/bin/benchmarks/advection_reaction_3D` directory:

* `advection_reaction_3D` -- MPI parallelism
* `advection_reaction_3D_mpicuda` -- MPI + CUDA parallelism
* `advection_reaction_3D_mpihip` -- MPI + HIP parallelism

On Summit, with the default environment
```
  Compiler: xl/16.1.1-5
  MPI: spectrum-mpi/10.3.1.2-20200121
  CUDA: cuda/10.1.243
```
an example `jsrun` command is
```
jsrun -n 2 -a 1 -c 1 -g 1 ./advection_reaction_3D_mpicuda
```

On Lassen, with the environment
```
  Compiler: gcc/8.3.1
  MPI: mvapich2/2021.05.28-cuda-11.1.1
  CUDA: cuda/11.1.1
```
an example `jsrun` command is
```
jsrun -n 2 -a 1 -c 1 -g 1 ./advection_reaction_3D_mpicuda
```
