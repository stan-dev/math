.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNNonlinSol:

###########################
Nonlinear Algebraic Solvers
###########################

SUNDIALS time integration packages are written in terms of generic nonlinear
solver operations defined by the SUNNonlinSol API and implemented by a
particular SUNNonlinSol module of type ``SUNNonlinearSolver``.
Users can supply their own SUNNonlinSol module, or use one of the modules
provided with SUNDIALS. Depending on the package, nonlinear solver modules
can either target systems presented in a rootfinding (:math:`F(y) = 0`) or
fixed-point (:math:`G(y) = y`) formulation. For more information on the
formulation of the nonlinear system(s) in CVODE, see
:numref:`SUNNonlinSol.CVODE`.


The time integrators in SUNDIALS specify a default nonlinear solver module
and as such this chapter is intended for users that wish to use a non-default
nonlinear solver module or would like to provide their own nonlinear solver
implementation. Users interested in using a non-default solver module may skip
the description of the SUNNonlinSol API in section :numref:`SUNNonlinSol.API`
and proceeded to the subsequent sections in this chapter that describe the
SUNNonlinSol modules provided with SUNDIALS.


For users interested in providing their own SUNNonlinSol module, the
following section presents the SUNNonlinSol API and its implementation
beginning with the definition of SUNNonlinSol functions in the
sections :numref:`SUNNonlinSol.API.CoreFn`, :numref:`SUNNonlinSol.API.SetFn` and
:numref:`SUNNonlinSol.API.GetFn`. This is followed by the definition of
functions supplied to a nonlinear solver implementation in the section
:numref:`SUNNonlinSol.API.SUNSuppliedFn`.  The nonlinear solver return
codes are given in the section :numref:`SUNNonlinSol.API.ReturnCodes`. The
``SUNNonlinearSolver`` type and the generic SUNNonlinSol module are defined
in the section :numref:`SUNNonlinSol.API.Generic`. Finally, the section
:numref:`SUNNonlinSol.API.Custom` lists the requirements for supplying a custom
SUNNonlinSol module. Users wishing to supply their own SUNNonlinSol module
are encouraged to use the SUNNonlinSol implementations provided with
SUNDIALS as templates for supplying custom nonlinear solver modules.




.. toctree::
   :maxdepth: 1

   SUNNonlinSol_API_link.rst
   CVODE_interface
   SUNNonlinSol_links.rst
