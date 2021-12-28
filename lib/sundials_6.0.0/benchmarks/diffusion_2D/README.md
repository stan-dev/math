# Benchmark: 2D Diffusion

This benchmark problem implements a 2D diffusion equation using either MPI,
MPI + CUDA, or MPI + HIP parallelism. Note a GPU-aware MPI implementation is
required.

## Problem description

This code simulates the anisotropic 2D heat equation,
```
    u_t = kx u_xx + ky u_yy + b,
```
where `kx` and `ky` are the diffusion coefficients. The system is evolved for
`t` in `[0, tf]` and `(x,y) = X` in `[0, X_max]^2` with the initial condition
```
    u(0,X) = sin^2(pi x) sin^2(pi y),
```
and stationary boundary conditions
```
     u_t(t,0,y) = u_t(t,x_max,y) = u_t(t,x,0) = u_t(t,x,y_max) = 0.
```
The source term is given by
```
     b(t,X) = -2 pi sin^2(pi x) sin^2(pi y) sin(pi t) cos(pi t)
              - kx 2 pi^2 (cos^2(pi x) - sin^2(pi x)) sin^2(pi y) cos^2(pi t)
              - ky 2 pi^2 (cos^2(pi y) - sin^2(pi y)) sin^2(pi x) cos^2(pi t).
```
Under this setup, the problem has the analytical solution
```
     u(t,X) = sin^2(pi x) sin^2(pi y) cos^2(pi t).
```
Spatial derivatives are computed using second-order centered differences on a
uniform spatial grid. The problem can be evolved in time with ARKODE, CVODE, or
IDA. With ARKODE, an adaptive step diagonally implicit Runge-Kutta (DIRK) method
is applied. When using CVODE or IDA, adaptive order and step BDF methods are
used.

In all cases, the nonlinear system(s) in each time step are solved using an
inexact Newton method paired with a matrix-free PCG or GMRES linear solver and a
Jacobi preconditioner.

## Options

Several command line options are available to change the problem parameters
as well as the integrator and solver options. A summary of the options are
listed below.

| Option                               | Description                                                                              | Default |
|:-------------------------------------|:-----------------------------------------------------------------------------------------|:--------|
| `--help`                             | Print the command line options and description                                           | --      |
| Problem Configuration Options        |                                                                                          |         |
| `--npx <int>`                        | Number of MPI tasks in the x-direction (0 forces MPI to decide)                          | 0       |
| `--npy <int>`                        | Number of MPI tasks in the y-direction (0 forces MPI to decide)                          | 0       |
| `--nx <int>`                         | Number of mesh points in the x-direction                                                 | 32      |
| `--ny <int>`                         | Number of mesh points in the y-direction                                                 | 32      |
| `--ux <realtype>`                    | The domain upper bound in the x-direction `x_max`                                        | 1.0     |
| `--uy <realtype>`                    | The domain upper bound in the y-direction `y_max`                                        | 1.0     |
| `--kx <realtype>`                    | Diffusion coefficient in the x-direction `kx`                                            | 1.0     |
| `--ky <realtype>`                    | Diffusion coefficient in the y-direction `ky`                                            | 1.0     |
| `--tf <realtype>`                    | The final time `tf`                                                                      | 1.0     |
| `--noforcing`                        | Disable the forcing term                                                                 | Enabled |
| Output Options                       |                                                                                          |         |
| `--output <int>`                     | Output level: `0` no output, `1` output progress and stats, `2` write solution to disk   | 1       |
| `--nout <int>`                       | Number of output times                                                                   | 20      |
| Common Integrator and Solver Options |                                                                                          |         |
| `--rtol <realtype>`                  | Relative tolerance                                                                       | 1e-5    |
| `--atol <realtype>`                  | Absolute tolerance                                                                       | 1e-10   |
| `--maxsteps <int>`                   | Max number of steps between outputs (0 uses the integrator default)                      | 0       |
| `--onstep <int>`                     | Number of steps to run using `ONE_STEP` mode for debugging (0 uses `NORMAL` mode)        | 0       |
| `--gmres`                            | Use GMRES rather than PCG                                                                | PCG     |
| `--lsinfo`                           | Output linear solver diagnostics                                                         | Off     |
| `--liniters <int>`                   | Number of linear iterations                                                              | 20      |
| `--epslin <realtype>`                | Linear solve tolerance factor (0 uses the integrator default)                            | 0       |
| `--msbp <int>`                       | The linear solver setup frequency (CVODE and ARKODE only, 0 uses the integrator default) | 0       |
| Additional ARKODE Options            |                                                                                          |         |
| `--order <int>`                      | Methods order                                                                            | 3       |
| `--controller <int>`                 | Error controller option                                                                  | 0       |
| `--nonlinear`                        | Treat the problem as nonlinearly implicit                                                | Linear  |
| `--diagnostics`                      | Output integrator diagnostics                                                            | Off     |

## Running

To build the benchmark executables SUNDIALS should be configured with ARKODE,
CVODE, or IDA enabled and with MPI support on. Additionally, either CUDA or HIP
support must be on to build executables utilizing NVIDIA or AMD GPUs. See the
installation guide for more details on configuring, building, and installing
SUNDIALS.

Based on the configuration, executables for each integrator and backend option
are built and installed in the `<install prefix>/bin/benchmarks/diffusion_2D`
directory. The executables follow the naming convention
`<package>_diffusion_2D_<parallelism>` where `<package>` is `arkode`,
`cvode`, or `ida` and `<parallelism>` is `mpi` for MPI only parallelism,
`mpicuda` for MPI + CUDA, and `mpihip` for MPI + HIP.

On Summit, with the default environment
```
  Compiler: xl/16.1.1-5
  MPI: spectrum-mpi/10.3.1.2-20200121
  CUDA: cuda/10.1.243
```
an example `jsrun` command using CUDA-aware MPI is
```
jsrun --smpiargs="-gpu" -n 2 -a 1 -c 1 -g 1 ./cvode_diffusion_2D_mpicuda
```

On Lassen, with the environment
```
  Compiler: gcc/8.3.1
  MPI: mvapich2/2021.05.28-cuda-11.1.1
  CUDA: cuda/11.1.1
```
an example `jsrun` command using CUDA-aware MPI is
```
jsrun -n 2 -a 1 -c 1 -g 1 ./cvode_diffusion_2D_mpicuda
```
