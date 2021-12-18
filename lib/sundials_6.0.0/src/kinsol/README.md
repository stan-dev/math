# KINSOL
### Version 6.0.0 (Dec 2021)

**Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, 
  and Carol S. Woodward, Center for Applied Scientific Computing, LLNL**

**Daniel R. Reynolds, Department of Mathematics, Southern Methodist University**


KINSOL is a package for the solution for nonlinear algebraic systems
```
F(u) = 0.
```
Nonlinear solver methods available include Newton-Krylov, Picard, and
fixed-point. Both Picard and fixed point can be accelerated with Anderson
acceleration.

KINSOL is part of the SUNDIALS Suite of Nonlinear and Differential/Algebraic
equation Solvers which consists of ARKode, CVODE, CVODES, IDA, IDAS and KINSOL.
It is written in ANSI standard C, but is based on the previous Fortran package
NKSOL, written by Peter Brown and Youcef Saad. KINSOL can be used in a variety
of computing environments including serial, shared memory, distributed memory,
and accelerator-based (e.g., GPU) systems. This flexibility is obtained from a
modular design that leverages the shared vector, matrix, and linear solver APIs
used across SUNDIALS packages.

For use with Fortran applications, a set of Fortran/C interface routines, called
FKINSOL, is also supplied.  These are written in C, but assume that the user
calling program and all user-supplied routines are in Fortran.

## Documentation

See the [KINSOL User Guide](/doc/kinsol/kin_guide.pdf) and
[KINSOL Examples](/doc/kinsol/kin_examples.pdf) document for more information
about IDA usage and the provided example programs respectively.

## Installation

For installation instructions see the [INSTALL_GUIDE](/INSTALL_GUIDE.pdf)
or the "Installation Procedure" chapter in the KINSOL User Guide.

## Release History

Information on recent changes to KINSOL can be found in the "Introduction"
chapter of the KINSOL User Guide and a complete release history is available in
the "SUNDIALS Release History" appendix of the KINSOL User Guide.

## References

* A. C. Hindmarsh, R. Serban, C. J. Balos, D. J. Gardner, 
  D. R. Reynolds and C. S. Woodward,
  "User Documentation for KINSOL v6.0.0," LLNL technical report
  UCRL-SM-208116, Dec 2021.

* A. M. Collier and R. Serban, "Example Programs for KINSOL v6.0.0,"
  LLNL technical report UCRL-SM-208114, Dec 2021.

* A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban,
  D. E. Shumaker, and C. S. Woodward, "SUNDIALS, Suite of Nonlinear and
  Differential/Algebraic Equation Solvers," ACM Trans. Math. Softw.,
  31(3), pp. 363-396, 2005.

* Peter N. Brown and Youcef Saad, "Hybrid Krylov Methods for
  Nonlinear Systems of Equations," SIAM J. Sci. Stat. Comput.,
  Vol 11, no 3, pp. 450-481, May 1990.

* A. G. Taylor and A. C. Hindmarsh, "User Documentation for KINSOL,
  A Nonlinear Solver for Sequential and Parallel Computers," LLNL
  technical report UCRL-ID-131185, July 1998.
