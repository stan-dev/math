# SUNDIALS: SUite of Nonlinear and DIfferential/ALgebraic equation Solvers #
### Version 5.2.0 (Mar 2020) ###

**Center for Applied Scientific Computing, Lawrence Livermore National Laboratory**

SUNDIALS is a family of software packages implemented with the goal of
providing robust time integrators and nonlinear solvers that can easily be
incorporated into existing simulation codes. The primary design goals are to
require minimal information from the user, allow users to easily supply their
own data structures underneath the packages, and allow for easy incorporation
of user-supplied linear solvers and preconditioners. The various packages share
many subordinate modules and are organized as a family with a directory
structure that exploits sharing common functionality.

The SUNDIALS suite consists of the following packages:

* ARKODE - for integration of stiff, nonstiff, and multirate ordinary
differential equation systems (ODEs) of the form

  ``` M y' = f1(t,y) + f2(t,y), y(t0) = y0 ```

* CVODE - for integration of stiff and nonstiff ordinary differential equation
systems (ODEs) of the form

  ``` y' = f(t,y), y(t0) = y0 ```

* CVODES - for integration and sensitivity analysis (forward and adjoint) of
ordinary differential equation systems (ODEs) of the form

  ``` y' = f(t,y,p), y(t0) = y0(p) ```

* IDA - for integration of differential-algebraic equation systems (DAEs) of
the form

  ``` F(t,y,y') = 0, y(t0) = y0, y'(t0) = y0' ```

* IDAS - for integration and sensitivity analysis (forward and adjoint) of
differential-algebraic equation systems (DAEs) of the form

  ``` F(t,y,y',p) = 0, y(t0) = y0(p), y'(t0) = y0'(p) ```

* KINSOL - for solution of nonlinear algebraic systems of the form

  ``` F(u) = 0 ```

## Installation ##

For installation directions see the [INSTALL_GUIDE](./INSTALL_GUIDE.pdf) or
the installation chapter in any of the package user guides.

Warning to users who receive more than one of the individual packages at
different times: Mixing old and new versions of SUNDIALS may fail. To avoid
such failures, obtain all desired package at the same time.

## Support ##

Full user guides for SUNDIALS packages are provided in the [doc](./doc)
directory along with documentation for example programs.

A list of Frequently Asked Questions on build and installation procedures as
well as common usage issues is available on the SUNDIALS [FAQ](https://computing.llnl.gov/projects/sundials/faq).
For dealing with systems with unphysical solutions or discontinuities see the
SUNDIALS [usage notes](https://computing.llnl.gov/projects/sundials/usage-notes).

If you have a question not covered in the FAQ or usage notes, please submit
your question to the SUNDIALS [mailing list](https://computing.llnl.gov/projects/sundials/mailing-list).

## Contributing ##

Bug fixes or minor changes are preferred via a pull request to the
[SUNDIALS GitHub repository](https://github.com/LLNL/sundials). For more
information on contributing see the [CONTRIBUTING](./CONTRIBUTING.md) file.

## Release History ##

Information on recent changes to SUNDIALS can be found in the "Introduction"
chapter of each package's user guide and a complete release history is available
in the "SUNDIALS Release History" appendix of each user guide.

## Authors ##

The SUNDIALS Team: Carol S. Woodward, Daniel R. Reynolds, Alan C. Hindmarsh,
David J. Gardner, and Cody J. Balos. We thank Radu Serban for significant
and critical past contributions.

We also acknowledge past contributions of Scott D. Cohen, Peter N. Brown,
George Byrne, Allan G. Taylor, Steven L. Lee, Keith E. Grant, Aaron Collier,
Lawrence E. Banks, Steve Smith, Cosmin Petra, Slaven Peles, John Loffeld, Dan Shumaker,
Ulrike Yang, James Almgren-Bell, Shelby L. Lockhart, Hilari C. Tiedeman, Ting Yan,
Jean M. Sexton, and Chris White.

### Citing SUNDIALS ###

We ask users of SUNDIALS to cite the following paper in any publications
reporting work done with SUNDIALS:

* Alan C. Hindmarsh, Peter N. Brown, Keith E. Grant, Steven L. Lee, Radu
Serban, Dan E. Shumaker, and Carol S. Woodward. 2005. SUNDIALS: Suite of
nonlinear and differential/algebraic equation solvers. ACM Trans. Math. Softw.
31, 3 (September 2005), 363-396. DOI=http://dx.doi.org/10.1145/1089014.1089020

## License ##

SUNDIALS is released under the BSD 3-clause license. See the [LICENSE](./LICENSE)
and [NOTICE](./NOTICE) files for details. All new contributions must be made
under the BSD 3-clause license.

**Please Note** If you are using SUNDIALS with any third party libraries linked
in (e.g., LAPACK, KLU, SuperLU_MT, PETSc, or *hypre*), be sure to review the
respective license of the package as that license may have more restrictive
terms than the SUNDIALS license.

```
SPDX-License-Identifier: BSD-3-Clause

LLNL-CODE-667205  (ARKODE)
UCRL-CODE-155951  (CVODE)
UCRL-CODE-155950  (CVODES)
UCRL-CODE-155952  (IDA)
UCRL-CODE-237203  (IDAS)
LLNL-CODE-665877  (KINSOL)
```
