# SUNDIALS: SUite of Nonlinear and DIfferential/ALgebraic equation Solvers #
### Version 6.1.1 (Feb 2022) ###

**Center for Applied Scientific Computing, Lawrence Livermore National Laboratory**

SUNDIALS is a family of software packages providing robust and efficient time
integrators and nonlinear solvers that can easily be incorporated into existing
simulation codes. The packages are designed to require minimal information from
the user, allow users to supply their own data structures underneath the
packages, and enable interfacing with user-supplied or third-party algebraic
solvers and preconditioners.

The SUNDIALS suite consists of the following packages for ordinary differential
equation (ODE) systems, differential-algebraic equation (DAE) systems, and
nonlinear algebraic systems:

* ARKODE - for integrating stiff, nonstiff, and multirate ODEs of the form

  `M(t) y' = f1(t,y) + f2(t,y), y(t0) = y0`

* CVODE - for integrating stiff and nonstiff ODEs of the form

  `y' = f(t,y), y(t0) = y0`

* CVODES - for integrating and sensitivity analysis (forward and adjoint) of
  ODEs of the form

  `y' = f(t,y,p), y(t0) = y0(p)`

* IDA - for integrating DAEs of the form

  `F(t,y,y') = 0, y(t0) = y0, y'(t0) = y0'`

* IDAS - for integrating and sensitivity analysis (forward and adjoint) of DAEs
  of the form

  `F(t,y,y',p) = 0, y(t0) = y0(p), y'(t0) = y0'(p)`

* KINSOL - for solving nonlinear algebraic systems of the form

  `F(u) = 0` or `G(u) = u`

## Installation ##

For installation directions see the [online documentation](https://sundials.readthedocs.io/en/latest/Install_link.html),
the [INSTALL_GUIDE](./INSTALL_GUIDE.pdf), or the installation chapter in any of
the package user guides.

Warning to users who receive more than one of the individual packages at
different times: Mixing old and new versions of SUNDIALS may fail. To avoid
such failures, obtain all desired package at the same time.

## Support ##

Full user guides for all of the SUNDIALS packages are available [online](https://sundials.readthedocs.io)
and in the [doc](./doc) directory. Additionally, the [doc](./doc) directory
contains documentation for the package example programs.

For information on recent changes to SUNDIALS see the [CHANGELOG](./CHANGELOG.md)
or the introduction chapter of any package user guide.

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

## Citing ##

See the [online documentation](https://sundials.readthedocs.io/en/latest/index.html#citing)
or [CITATIONS](./CITATIONS.md) file for information on how to cite SUNDIALS in
any publications reporting work done using SUNDIALS packages.

## Authors ##

The SUNDIALS library has been developed over many years by a number of
contributors. The current SUNDIALS team consists of Cody J. Balos,
David J. Gardner, Alan C. Hindmarsh, Daniel R. Reynolds, and Carol S. Woodward.
We thank Radu Serban for significant and critical past contributions.

Other contributors to SUNDIALS include: James Almgren-Bell, Lawrence E. Banks,
Peter N. Brown, George Byrne, Rujeko Chinomona, Scott D. Cohen, Aaron Collier,
Keith E. Grant, Steven L. Lee, Shelby L. Lockhart, John Loffeld, Daniel McGreer,
Slaven Peles, Cosmin Petra, H. Hunter Schwartz, Jean M. Sexton,
Dan Shumaker, Steve G. Smith, Allan G. Taylor, Hilari C. Tiedeman, Chris White,
Ting Yan, and Ulrike M. Yang.

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
