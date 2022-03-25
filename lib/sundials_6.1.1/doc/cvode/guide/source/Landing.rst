.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. only:: html

   This is the documentation for CVODE, a component of the `SUNDIALS
   <https://computing.llnl.gov/casc/sundials/main.html>`_ suite of
   nonlinear and differential/algebraic equation solvers.

   This work was performed under the auspices of the U.S. Department of Energy by Lawrence
   Livermore National Laboratory under Contract DE-AC52-07NA27344.

   **Authors:** *Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, Daniel R. Reynolds, and Carol S. Woodward.*

   **Contributors:** The SUNDIALS library has been developed over many years by
   a number of contributors. The current SUNDIALS team consists of Cody J.
   Balos, David J. Gardner, Alan C. Hindmarsh, Daniel R. Reynolds, and Carol S.
   Woodward. We thank Radu Serban for significant and critical past
   contributions.

   Other contributors to SUNDIALS include: James Almgren-Bell, Lawrence E.
   Banks, Peter N. Brown, George Byrne, Rujeko Chinomona, Scott D. Cohen, Aaron
   Collier, Keith E. Grant, Steven L. Lee, Shelby L. Lockhart, John Loffeld,
   Daniel McGreer, Slaven Peles, Cosmin Petra, H. Hunter Schwartz, Jean M.
   Sexton, Dan Shumaker, Steve G. Smith, Allan G. Taylor, Hilari C. Tiedeman,
   Chris White, Ting Yan, and Ulrike M. Yang.

   .. ifconfig:: package_name != 'super'

      **Citing**

      .. include:: ../../../shared/cite_sundials.rst

      The CVODE documentation can be cited:

      .. parsed-literal::

         @Misc{cvodeDocumentation,
            author = {Alan C. Hindmarsh and Radu Serban and Cody J. Balos and David J. Gardner and Daniel R. Reynolds and Carol S. Woodward},
            title  = {User Documentation for CVODE},
            year   = {|YEAR|}
            notes  = {|CVODE_VERSION|}
         }
