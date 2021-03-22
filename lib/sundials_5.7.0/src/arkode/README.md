# ARKode
### Version 4.7.0 (Jan 2021)

**Daniel R. Reynolds,
  Department of Mathematics, SMU**

**David J. Gardner, Carol S. Woodward, and Cody J. Balos,  
  Center for Applied Scientific Computing, LLNL**

ARKode is a package for the solution of stiff, nonstiff, and multirate ordinary
differential equation (ODE) systems (initial value problems) given in linearly
implicit the form
```
M y' = f1(t,y) + f2(t,y), y(t0) = y0.
```
The integration methods implemented in ARKode include explicit and implicit
Runge-Kutta methods, implicit-explicit (IMEX) additive Runge-Kutta methods, and
multirate infatesemial step (MIS) methods.

ARKode is part of a the SUNDIALS Suite of Nonlinear and Differential/Algebraic
equation Solvers which consists of ARKode, CVODE, CVODES, IDA, IDAS, and KINSOL.
It is written in ANSI standard C and can be used in a variety of computing
environments including serial, shared memory, distributed memory, and
accelerator-based (e.g., GPU) systems. This flexibility is obtained from a
modular design that leverages the shared vector, matrix, linear solver, and
nonlinear solver APIs used across SUNDIALS packages.

For use with Fortran applications, a set of Fortran/C interface routines, called
FARKODE, is also supplied. These are written in C, but assume that the user
calling program and all user-supplied routines are in Fortran.

## Documentation

See the [ARKode User Guide](/doc/arkode/ark_guide.pdf) and
[ARKode Examples](/doc/arkode/ark_examples.pdf) document for more information
about ARKode usage and the provided example programs respectively.

## Installation

For installation instructions see the [INSTALL_GUIDE](/INSTALL_GUIDE.pdf)
or "Installation Procedure" chapter in the ARKode User Guide.

## Release History

Information on recent changes to ARKode can be found in the "Introduction"
chapter of the ARKode User Guide and a complete release history is available in
the "SUNDIALS Release History" appendix of the ARKode User Guide.

## References

* D. R. Reynolds, D. J. Gardner, C. S. Woodward, and C. J. Balos,
  "User Documentation for ARKode v4.7.0," LLNL technical report
  LLNL-SM-668082, Jan 2021.

* D. R. Reynolds, "Example Programs for ARKode v4.7.0," Technical Report,
  Southern Methodist University Center for Scientific Computation, Jan 2021.
