# CVODES
### Version 6.0.0 (Dec 2021)

**Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, 
  and Carol S. Woodward, Center for Applied Scientific Computing, LLNL**

**Daniel R. Reynolds, Department of Mathematics, Southern Methodist University**


CVODES is a package for the solution of stiff and nonstiff ordinary differential
equation (ODE) systems (initial value problem) given in explicit form
```
dy/dt = f(t,y,p), y(t0) = y0(p)
```
with sensitivity analysis capabilities (both forward and adjoint modes). CVODES
provides a choice of two variable-order, variable-coefficient multistep methods,
Adams-Moulton methods for non-stiff problems or BDF (Backward Differentiation
Formula) methods in fixed-leading-coefficient form for stiff problems.

CVODES is part of the SUNDIALS Suite of Nonlinear and Differential/Algebraic
equation Solvers which consists of ARKode, CVODE, CVODES, IDA, IDAS and KINSOL.
It is written in ANSI standard C and can be used in a variety of computing
environments including serial, shared memory, distributed memory, and
accelerator-based (e.g., GPU) systems. This flexibility is obtained from a
modular design that leverages the shared vector, matrix, linear solver, and
nonlinear solver APIs used across SUNDIALS packages.

## Documentation

See the [CVODES User Guide](/doc/cvodes/cvs_guide.pdf) and
[CVODES Examples](/doc/cvodes/cvs_examples.pdf) document for more information
about CVODES usage and the provided example programs respectively.

## Installation

For installation instructions see the [INSTALL_GUIDE](/INSTALL_GUIDE.pdf)
or the "Installation Procedure" chapter in the CVODES User Guide.

## Release History

Information on recent changes to CVODES can be found in the "Introduction"
chapter of the CVODES User Guide and a complete release history is available in
the "SUNDIALS Release History" appendix of the CVODES User Guide.

## References

* A. C. Hindmarsh, R. Serban, C. J. Balos, D. J. Gardner, D. R. Reynolds
  and C. S. Woodward, "User Documentation for CVODES v6.0.0,"
  LLNL technical report UCRL-SM-208111, Dec 2021.

* A. C. Hindmarsh and R. Serban, "Example Programs for CVODES v6.0.0,"
  LLNL technical report UCRL-SM-208115, Dec 2021.

* R. Serban and A. C. Hindmarsh, "CVODES: the Sensitivity-Enabled ODE
  solver in SUNDIALS," Proceedings of IDETC/CIE 2005, Sept. 2005,
  Long Beach, CA.

* A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban,
  D. E. Shumaker, and C. S. Woodward, "SUNDIALS, Suite of Nonlinear and
  Differential/Algebraic Equation Solvers," ACM Trans. Math. Softw.,
  31(3), pp. 363-396, 2005.
