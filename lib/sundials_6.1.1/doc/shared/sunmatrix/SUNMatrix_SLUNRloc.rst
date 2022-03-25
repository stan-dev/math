..
   Programmer(s): Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMatrix.SLUNRloc:

The SUNMATRIX_SLUNRLOC Module
======================================

The SUNMATRIX_SLUNRLOC module is an interface to the ``SuperMatrix``
structure provided by the SuperLU_DIST sparse matrix factorization and
solver library written by X. Sherry Li and collaborators
:cite:p:`SuperLUDIST_site,GDL:07,LD:03,SLUUG:99`.
It is designed to be used with the SuperLU_DIST ``SUNLinearSolver``
module discussed in :numref:`SUNLinSol.SuperLUDIST`. To this end, it
defines the ``content`` field of ``SUNMatrix`` to be the following
structure:

.. code-block:: c

   struct _SUNMatrixContent_SLUNRloc {
     booleantype   own_data;
     gridinfo_t    *grid;
     sunindextype  *row_to_proc;
     pdgsmv_comm_t *gsmv_comm;
     SuperMatrix   *A_super;
     SuperMatrix   *ACS_super;
   };

A more complete description of the this ``content`` field is given below:

* ``own_data`` -- a flag which indicates if the SUNMatrix is responsible for freeing
  ``A_super``

* ``grid`` -- pointer to the SuperLU_DIST structure that stores the 2D process grid

* ``row_to_proc`` -- a mapping between the rows in the matrix and the process it
  resides on; will be ``NULL`` until the ``SUNMatMatvecSetup`` routine is called

* ``gsmv_comm`` -- pointer to the SuperLU_DIST structure that stores the
  communication information needed for matrix-vector multiplication; will be
  ``NULL`` until the ``SUNMatMatvecSetup`` routine is called

* ``A_super`` -- pointer to the underlying SuperLU_DIST ``SuperMatrix`` with
  ``Stype = SLU_NR_loc``, ``Dtype = SLU_D``, ``Mtype = SLU_GE``; must have the
  full diagonal present to be used with ``SUNMatScaleAddI`` routine

* ``ACS_super`` -- a column-sorted version of the matrix needed to perform matrix-vector
  multiplication; will be ``NULL`` until the routine ``SUNMatMatvecSetup``
  routine is called


The header file to include when using this module is ``sunmatrix/sunmatrix_slunrloc.h``.
The installed module library to link to is ``libsundials_sunmatrixslunrloc`` *.lib*
where *.lib* is typically ``.so`` for shared libraries and ``.a`` for static libraries.


.. _SUNMatrix.SLUNRloc.Functions:

SUNMATRIX_SLUNRLOC Functions
----------------------------------

The SUNMATRIX_SLUNRLOC module provides the following user-callable routines:


.. c:function:: SUNMatrix SUNMatrix_SLUNRloc(SuperMatrix *Asuper, gridinfo_t *grid, SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNMATRIX_SLUNRLOC
   object. Its arguments are a fully-allocated SuperLU_DIST ``SuperMatrix`` with
   ``Stype = SLU_NR_loc, Dtype = SLU_D, Mtype = SLU_GE`` and an initialized SuperLU_DIST
   2D process grid structure. It returns a ``SUNMatrix`` object if ``Asuper`` is
   compatible else it returns ``NULL``.


.. c:function:: void SUNMatrix_SLUNRloc_Print(SUNMatrix A, FILE *fp)

   This function prints the underlying ``SuperMatrix`` content. It is useful for
   debugging. Its arguments are the ``SUNMatrix`` object and a ``FILE`` pointer
   to print to. It returns void.


.. c:function:: SuperMatrix* SUNMatrix_SLUNRloc_SuperMatrix(SUNMatrix A)

   This function returns the underlying ``SuperMatrix`` of ``A``. Its
   only argument is the ``SUNMatrix`` object to access.


.. c:function:: gridinfo_t* SUNMatrix_SLUNRloc_ProcessGrid(SUNMatrix A)

   This function returns the SuperLU_DIST 2D process grid associated with
   ``A``. Its only argument is the ``SUNMatrix`` object to access.


.. c:function:: booleantype SUNMatrix_SLUNRloc_OwnData(SUNMatrix A)

   This function returns true if the ``SUNMatrix`` object is responsible
   for freeing the underlying ``SuperMatrix``, otherwise it returns false.
   Its only argument is the ``SUNMatrix`` object to access.


The SUNMATRIX_SLUNRLOC module also defines implementations of all generic
``SUNMatrix`` operations listed in :numref:`SUNMatrix.ops`:

* ``SUNMatGetID_SLUNRloc`` -- returns ``SUNMATRIX_SLUNRLOC``

* ``SUNMatClone_SLUNRloc``

* ``SUNMatDestroy_SLUNRloc``

* ``SUNMatSpace_SLUNRloc`` -- this only returns information for the storage within
  the matrix interface, i.e. storage for ``row_to_proc``

* ``SUNMatZero_SLUNRloc``

* ``SUNMatCopy_SLUNRloc``

* ``SUNMatScaleAdd_SLUNRloc`` -- performs :math:`A = cA + B`, where :math:`A` and :math:`B`
  must have the same sparsity pattern

* ``SUNMatScaleAddI_SLUNRloc`` -- performs :math:`A = cA + I`, where the diagonal of :math:`A`
  must be present

* ``SUNMatMatvecSetup_SLUNRloc`` -- initializes the SuperLU_DIST parallel communication
  structures needed to perform a matrix-vector product; only needs to be called before
  the first call to :c:func:`SUNMatMatvec` or if the matrix changed since the last setup

* ``SUNMatMatvec_SLUNRloc``
