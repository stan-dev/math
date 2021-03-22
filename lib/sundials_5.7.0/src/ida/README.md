# IDA
### Version 5.7.0 (Jan 2021)

**Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, 
  and Carol S. Woodward, Center for Applied Scientific Computing, LLNL**

**Daniel R. Reynolds, Department of Mathematics, Southern Methodist University**


IDA is a package for the solution of differential-algebraic equation (DAE)
systems
```
F(t,y,y') = 0, y(t0) = y0, y'(t0) = y0'.
```
The integration methods used in IDA are variable-order, variable-coefficient BDF
(Backward Differentiation Formula) methods in fixed-leading-coefficient form.

IDA is part of the SUNDIALS Suite of Nonlinear and Differential/Algebraic
equation Solvers which consists of ARKode, CVODE, CVODES, IDA, IDAS and KINSOL.
It is written in ANSI standard C, but is derived from the package DASPK, which
is written in FORTRAN. IDA can be used in a variety of computing environments
including serial, shared memory, distributed memory, and accelerator-based
(e.g., GPU) systems. This flexibility is obtained from a modular design that
leverages the shared vector, matrix, linear solver, and nonlinear solver APIs
used across SUNDIALS packages.

For use with Fortran applications, a set of Fortran/C interface routines, called
FIDA, is also supplied. These are written in C, but assume that the user calling
program and all user-supplied routines are in Fortran.

## Documentation

See the [IDA User Guide](/doc/ida/ida_guide.pdf) and
[IDA Examples](/doc/ida/ida_examples.pdf) document for more information
about IDA usage and the provided example programs respectively.

## Installation

For installation instructions see the [INSTALL_GUIDE](/INSTALL_GUIDE.pdf)
or the "Installation Procedure" chapter in the IDA User Guide.

## Release History

Information on recent changes to IDA can be found in the "Introduction"
chapter of the IDA User Guide and a complete release history is available in
the "SUNDIALS Release History" appendix of the IDA User Guide.

## References

* A. C. Hindmarsh, R. Serban, C. J. Balos, D. J. Gardner, D. R. Reynolds
  and C. S. Woodward, "User Documentation for IDA v5.7.0,"
  LLNL technical report UCRL-SM-208112, Jan 2021.

* A. C. Hindmarsh, R. Serban, and A. Collier, "Example Programs for IDA v5.7.0,"
  LLNL technical report UCRL-SM-208113, Jan 2021.

* A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban,
  D. E. Shumaker, and C. S. Woodward, "SUNDIALS, Suite of Nonlinear and
  Differential/Algebraic Equation Solvers," ACM Trans. Math. Softw.,
  31(3), pp. 363-396, 2005.

* P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov Methods
  in the Solution of Large-Scale Differential-Algebraic Systems,
  SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488.

* P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent Initial
  Condition Calculation for Differential-Algebraic Systems,
  SIAM J. Sci. Comp., 19 (1998), pp. 1495-1512.
