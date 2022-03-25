.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS.Fortran:

SUNDIALS Fortran Interface
==========================

SUNDIALS provides modern, Fortran 2003 based, interfaces as Fortran modules to
most of the C API including:

- All of the time-stepping modules in ARKODE:

  * The ``farkode_arkstep_mod``, ``farkode_erkstep_mod``, and
    ``farkode_mristep_mod`` modules provide interfaces to the ARKStep, ERKStep,
    and MRIStep integrators respectively.

  * The ``farkode_mod`` module interfaces to the components of ARKODE which are
    shared by the time-stepping modules.

- CVODE via the ``fcvode_mod`` module.

- CVODES via the ``fcvodes_mod`` module.

- IDA via the ``fida_mod`` module.

- IDAS via the ``fidas_mod`` module.

- KINSOL via the ``fkinsol_mod`` module.

.. ifconfig:: package_name == 'kinsol'

   Additionally, all of the SUNDIALS base classes (:c:type:`N_Vector`, :c:type:`SUNMatrix`, and
   :c:type:`SUNLinearSolver`) include Fortran interface modules.  A complete list of class
   implementations with Fortran 2003 interface modules is given in :numref:`SUNDIALS.Fortran.Table`.

.. ifconfig:: package_name != 'kinsol'

   Additionally, all of the SUNDIALS base classes (:c:type:`N_Vector`, :c:type:`SUNMatrix`,
   :c:type:`SUNLinearSolver`, and :c:type:`SUNNonlinearSolver`) include Fortran interface modules.
   A complete list of class implementations with Fortran 2003 interface modules is given in
   :numref:`SUNDIALS.Fortran.Table`.

An interface module can be accessed with the ``use`` statement, e.g.

.. code-block:: fortran

   use fcvode_mod
   use fnvector_openmp_mod

and by linking to the Fortran 2003 library in addition to the C library, e.g.
``libsundials_fnvecpenmp_mod.<so|a>``, ``libsundials_nvecopenmp.<so|a>``,
``libsundials_fcvode_mod.<so|a>`` and ``libsundials_cvode.<so|a>``.

The Fortran 2003 interfaces leverage the ``iso_c_binding`` module and the
``bind(C)`` attribute to closely follow the SUNDIALS C API (modulo language
differences). The SUNDIALS classes, e.g. :c:type:`N_Vector`, are interfaced as
Fortran derived types, and function signatures are matched but with an ``F``
prepending the name, e.g. ``FN_VConst`` instead of :c:func:`N_VConst` or
``FCVodeCreate`` instead of ``CVodeCreate``. Constants are named exactly as they
are in the C API.  Accordingly, using SUNDIALS via the Fortran 2003 interfaces
looks just like using it in C. Some caveats stemming from the language
differences are discussed in :numref:`SUNDIALS.Fortran.Differences`. A
discussion on the topic of equivalent data types in C and Fortran 2003 is
presented in :numref:`SUNDIALS.Fortran.DataTypes`.

.. ifconfig:: package_name == 'kinsol'

   Further information on the Fortran 2003 interfaces specific to the
   :c:type:`N_Vector`, :c:type:`SUNMatrix`, and :c:type:`SUNLinearSolver`
   classes is given alongside the C documentation (:numref:`NVectors`,
   :numref:`SUNMatrix`, and :numref:`SUNLinSol`, respectively). For details on
   where the Fortran 2003 module (``.mod``) files and libraries are installed
   see :numref:`Installation`.

.. ifconfig:: package_name != 'kinsol'

   Further information on the Fortran 2003 interfaces specific to the
   :c:type:`N_Vector`, :c:type:`SUNMatrix`, :c:type:`SUNLinearSolver`, and
   :c:type:`SUNNonlinearSolver` classes is given alongside the C documentation
   (:numref:`NVectors`, :numref:`SUNMatrix`, :numref:`SUNLinSol`, and
   :numref:`SUNNonlinSol` respectively). For details on where the Fortran 2003
   module (``.mod``) files and libraries are installed see
   :numref:`Installation`.

The Fortran 2003 interface modules were generated with SWIG Fortran
:cite:p:`Swig-Fortran`, a fork of SWIG. Users who are interested in the SWIG
code used in the generation process should contact the SUNDIALS development
team.

.. _SUNDIALS.Fortran.Table:

.. table:: List of SUNDIALS Fortran 2003 interface modules

   =======================  ====================================
   **Class/Module**          **Fortran 2003 Module Name**
   =======================  ====================================
   ARKODE                   ``farkode_mod``
   ARKODE::ARKSTEP          ``farkode_arkstep_mod``
   ARKODE::ERKSTEP          ``farkode_erkstep_mod``
   ARKODE::MRISTEP          ``farkode_mristep_mod``
   CVODE                    ``fcvode_mod``
   CVODES                   ``fcvodes_mod``
   IDA                      ``fida_mod``
   IDAS                     ``fidas_mod``
   KINSOL                   ``fkinsol_mod``
   NVECTOR                  ``fsundials_nvector_mod``
   NVECTOR_SERIAL           ``fnvector_serial_mod``
   NVECTOR_OPENMP           ``fnvector_openmp_mod``
   NVECTOR_PTHREADS         ``fnvector_pthreads_mod``
   NVECTOR_PARALLEL         ``fnvector_parallel_mod``
   NVECTOR_PARHYP           Not interfaced
   NVECTOR_PETSC            Not interfaced
   NVECTOR_CUDA             Not interfaced
   NVECTOR_RAJA             Not interfaced
   NVECTOR_SYCL             Not interfaced
   NVECTOR_MANVECTOR        ``fnvector_manyvector_mod``
   NVECTOR_MPIMANVECTOR     ``fnvector_mpimanyvector_mod``
   NVECTOR_MPIPLUSX         ``fnvector_mpiplusx_mod``
   SUNMATRIX                ``fsundials_matrix_mod``
   SUNMATRIX_BAND           ``fsunmatrix_band_mod``
   SUNMATRIX_DENSE          ``fsunmatrix_dense_mod``
   SUNMATRIX_MAGMADENSE     Not interfaced
   SUNMATRIX_ONEMKLDENSE    Not interfaced
   SUNMATRIX_SPARSE         ``fsunmatrix_sparse_mod``
   SUNLINSOL                ``fsundials_linearsolver_mod``
   SUNLINSOL_BAND           ``fsunlinsol_band_mod``
   SUNLINSOL_DENSE          ``fsunlinsol_dense_mod``
   SUNLINSOL_LAPACKBAND     Not interfaced
   SUNLINSOL_LAPACKDENSE    Not interfaced
   SUNLINSOL_MAGMADENSE     Not interfaced
   SUNLINSOL_ONEMKLDENSE    Not interfaced
   SUNLINSOL_KLU            ``fsunlinsol_klu_mod``
   SUNLINSOL_SLUMT          Not interfaced
   SUNLINSOL_SLUDIST        Not interfaced
   SUNLINSOL_SPGMR          ``fsunlinsol_spgmr_mod``
   SUNLINSOL_SPFGMR         ``fsunlinsol_spfgmr_mod``
   SUNLINSOL_SPBCGS         ``fsunlinsol_spbcgs_mod``
   SUNLINSOL_SPTFQMR        ``fsunlinsol_sptfqmr_mod``
   SUNLINSOL_PCG            ``fsunlinsol_pcg_mof``
   SUNNONLINSOL             ``fsundials_nonlinearsolver_mod``
   SUNNONLINSOL_NEWTON      ``fsunnonlinsol_newton_mod``
   SUNNONLINSOL_FIXEDPOINT  ``fsunnonlinsol_fixedpoint_mod``
   SUNNONLINSOL_PETSCSNES   Not interfaced
   =======================  ====================================


.. _SUNDIALS.Fortran.DataTypes:

Data Types
----------

Generally, the Fortran 2003 type that is equivalent to the C type is what one
would expect. Primitive types map to the ``iso_c_binding`` type equivalent.
SUNDIALS classes map to a Fortran derived type. However, the handling of pointer
types is not always clear as they can depend on the parameter direction.
:numref:`SUNDIALS.Fortran.DataTypes.Table` presents a summary of the type
equivalencies with the parameter direction in mind.

.. warning::

   Currently, the Fortran 2003 interfaces are only compatible with SUNDIALS
   builds where the ``realtype`` is double-precision the ``sunindextype`` size
   is 64-bits.

.. _SUNDIALS.Fortran.DataTypes.Table:
.. table:: C/Fortran-2003 Equivalent Types

   +-------------------------+-------------------------------+-------------------------------------------+
   | **C Type**              | **Parameter Direction**       | **Fortran 2003 type**                     |
   +=========================+===============================+===========================================+
   |``double``               | in, inout, out, return        | ``real(c_double)``                        |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``int``                  | in, inout, out, return        | ``integer(c_int)``                        |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``long``                 | in, inout, out, return        | ``integer(c_long)``                       |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``booleantype``          | in, inout, out, return        | ``integer(c_int)``                        |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``realtype``             | in, inout, out, return        | ``real(c_double)``                        |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``sunindextype``         | in, inout, out, return        | ``integer(c_long)``                       |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``double*``              | in, inout, out                | ``real(c_double), dimension(*)``          |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``double*``              | return                        | ``real(c_double), pointer, dimension(:)`` |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``int*``                 | in, inout, out                | ``real(c_int), dimension(*)``             |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``int*``                 | return                        | ``real(c_int), pointer, dimension(:)``    |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``long*``                | in, inout, out                | ``real(c_long), dimension(*)``            |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``long*``                | return                        | ``real(c_long), pointer, dimension(:)``   |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``realtype*``            | in, inout, out                | ``real(c_double), dimension(*)``          |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``realtype*``            | return                        | ``real(c_double), pointer, dimension(:)`` |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``sunindextype*``        | in, inout, out                | ``real(c_long), dimension(*)``            |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``sunindextype*``        | return                        | ``real(c_long), pointer, dimension(:)``   |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``realtype[]``           | in, inout, out                | ``real(c_double), dimension(*)``          |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``sunindextype[]``       | in, inout, out                | ``integer(c_long), dimension(*)``         |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``N_Vector``             | in, inout, out                | ``type(N_Vector)``                        |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``N_Vector``             | return                        | ``type(N_Vector), pointer``               |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``SUNMatrix``            | in, inout, out                | ``type(SUNMatrix)``                       |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``SUNMatrix``            | return                        | ``type(SUNMatrix), pointer``              |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``SUNLinearSolver``      | in, inout, out                | ``type(SUNLinearSolver)``                 |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``SUNLinearSolver``      | return                        | ``type(SUNLinearSolver), pointer``        |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``SUNNonlinearSolver``   | in, inout, out                | ``type(SUNNonlinearSolver)``              |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``SUNNonlinearSolver``   | return                        | ``type(SUNNonlinearSolver), pointer``     |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``FILE*``                | in, inout, out, return        | ``type(c_ptr)``                           |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``void*``                | in, inout, out, return        | ``type(c_ptr)``                           |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``T**``                  | in, inout, out, return        | ``type(c_ptr)``                           |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``T***``                 | in, inout, out, return        | ``type(c_ptr)``                           |
   +-------------------------+-------------------------------+-------------------------------------------+
   |``T****``                | in, inout, out, return        | ``type(c_ptr)``                           |
   +-------------------------+-------------------------------+-------------------------------------------+


.. _SUNDIALS.Fortran.Differences:

Notable Fortran/C usage differences
-----------------------------------

While the Fortran 2003 interface to SUNDIALS closely follows the C API, some
differences are inevitable due to the differences between Fortran and C.  In
this section, we note the most critical differences. Additionally,
:numref:`SUNDIALS.Fortran.DataTypes` discusses equivalencies of data types
in the two languages.


.. _SUNDIALS.Fortran.Differences.CreatingObjects:

Creating generic SUNDIALS objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the C API a SUNDIALS class, such as an :c:type:`N_Vector`, is actually a pointer to
an underlying C struct. However, in the Fortran 2003 interface, the derived type
is bound to the C struct, not the pointer to the struct. For example,
``type(N_Vector)`` is bound to the C struct ``_generic_N_Vector`` not the
``N_Vector`` type. The consequence of this is that creating and declaring SUNDIALS
objects in Fortran is nuanced. This is illustrated in the code snippets below:

C code:

.. sourcecode:: c

   N_Vector x;
   x = N_VNew_Serial(N, sunctx);

Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x
   x => FN_VNew_Serial(N, sunctx)

Note that in the Fortran declaration, the vector is a ``type(N_Vector),
pointer``, and that the pointer assignment operator is then used.


.. _SUNDIALS.Fortran.Differences.ArraysAndPointers:

Arrays and pointers
^^^^^^^^^^^^^^^^^^^

Unlike in the C API, in the Fortran 2003 interface, arrays and pointers are
treated differently when they are return values versus arguments to a function.
Additionally, pointers which are meant to be out parameters, not arrays, in the
C API must still be declared as a rank-1 array in Fortran.  The reason for this
is partially due to the Fortran 2003 standard for C bindings, and partially due
to the tool used to generate the interfaces. Regardless, the code snippets below
illustrate the differences.

C code:

.. sourcecode:: c

   N_Vector x;
   realtype* xdata;
   long int leniw, lenrw;

   /* create a new serial vector */
   x = N_VNew_Serial(N, sunctx);

   /* capturing a returned array/pointer */
   xdata = N_VGetArrayPointer(x)

   /* passing array/pointer to a function */
   N_VSetArrayPointer(xdata, x)

   /* pointers that are out-parameters */
   N_VSpace(x, &leniw, &lenrw);


Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x
   real(c_double), pointer :: xdataptr(:)
   real(c_double)          :: xdata(N)
   integer(c_long)         :: leniw(1), lenrw(1)

   ! create a new serial vector
   x => FN_VNew_Serial(x, sunctx)

   ! capturing a returned array/pointer
   xdataptr => FN_VGetArrayPointer(x)

   ! passing array/pointer to a function
   call FN_VSetArrayPointer(xdata, x)

   ! pointers that are out-parameters
   call FN_VSpace(x, leniw, lenrw)


.. _SUNDIALS.Fortran.Differences.ProcedurePointers:

Passing procedure pointers and user data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since functions/subroutines passed to SUNDIALS will be called from within C
code, the Fortran procedure must have the attribute ``bind(C)``. Additionally,
when providing them as arguments to a Fortran 2003 interface routine, it is
required to convert a procedure's Fortran address to C with the Fortran
intrinsic ``c_funloc``.

Typically when passing user data to a SUNDIALS function, a user may simply cast
some custom data structure as a ``void*``. When using the Fortran 2003
interfaces, the same thing can be achieved. Note, the custom data structure
*does not* have to be ``bind(C)`` since it is never accessed on the C side.

C code:

.. sourcecode:: c

   MyUserData *udata;
   void *cvode_mem;

   ierr = CVodeSetUserData(cvode_mem, udata);

Fortran code:

.. sourcecode:: Fortran

   type(MyUserData) :: udata
   type(c_ptr)      :: arkode_mem

   ierr = FARKStepSetUserData(arkode_mem, c_loc(udata))

On the other hand, Fortran users may instead choose to store problem-specific
data, e.g.  problem parameters, within modules, and thus do not need the
SUNDIALS-provided ``user_data`` pointers to pass such data back to user-supplied
functions. These users should supply the ``c_null_ptr`` input for ``user_data``
arguments to the relevant SUNDIALS functions.

.. _SUNDIALS.Fortran.Differences.OptionalParameters:

Passing ``NULL`` to optional parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the SUNDIALS C API some functions have optional parameters that a caller can
pass as ``NULL``. If the optional parameter is of a type that is equivalent to a
Fortran ``type(c_ptr)`` (see :numref:`SUNDIALS.Fortran.DataTypes`),
then a Fortran user can pass the intrinsic ``c_null_ptr``. However, if the
optional parameter is of a type that is not equivalent to ``type(c_ptr)``, then
a caller must provide a Fortran pointer that is dissociated. This is
demonstrated in the code example below.

C code:

.. sourcecode:: c

   SUNLinearSolver LS;
   N_Vector x, b;

   /* SUNLinSolSolve expects a SUNMatrix or NULL as the second parameter. */
   ierr = SUNLinSolSolve(LS, NULL, x, b);

Fortran code:

.. sourcecode:: Fortran

   type(SUNLinearSolver), pointer :: LS
   type(SUNMatrix), pointer       :: A
   type(N_Vector), pointer        :: x, b

   ! Disassociate A
   A => null()

   ! SUNLinSolSolve expects a type(SUNMatrix), pointer as the second parameter.
   ! Therefore, we cannot pass a c_null_ptr, rather we pass a disassociated A.
   ierr = FSUNLinSolSolve(LS, A, x, b)

.. _SUNDIALS.Fortran.Differences.NVectorArrays:

Working with ``N_Vector`` arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Arrays of :c:type:`N_Vector` objects are interfaced to Fortran 2003 as an opaque
``type(c_ptr)``.  As such, it is not possible to directly index an array of
:c:type:`N_Vector` objects returned by the ``N_Vector`` "VectorArray" operations, or
packages with sensitivity capabilities (CVODES and IDAS).  Instead, SUNDIALS
provides a utility function :f:func:`FN_VGetVecAtIndexVectorArray` that can be
called for accessing a vector in a vector array. The example below demonstrates
this:

C code:

.. sourcecode:: c

   N_Vector x;
   N_Vector* vecs;

   /* Create an array of N_Vectors */
   vecs = N_VCloneVectorArray(count, x);

   /* Fill each array with ones */
   for (int i = 0; i < count; ++i)
     N_VConst(vecs[i], 1.0);

Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x, xi
   type(c_ptr)             :: vecs

   ! Create an array of N_Vectors
   vecs = FN_VCloneVectorArray(count, x)

   ! Fill each array with ones
   do index = 0,count-1
     xi => FN_VGetVecAtIndexVectorArray(vecs, index)
     call FN_VConst(xi, 1.d0)
   enddo

SUNDIALS also provides the functions :c:func:`N_VSetVecAtIndexVectorArray` and
:c:func:`N_VNewVectorArray` for working with ``N_Vector`` arrays, that have
corresponding Fortran interfaces ``FN_VSetVecAtIndexVectorArray`` and
``FN_VNewVectorArray``, respectively. These functions are particularly
useful for users of the Fortran interface to the
:ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` or
:ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` when creating the
subvector array. Both of these functions along with
:c:func:`N_VGetVecAtIndexVectorArray` (wrapped as
``FN_VGetVecAtIndexVectorArray``) are further described in
:numref:`NVectors.Description.utilities`.

.. _SUNDIALS.Fortran.Differences.FilePointers:

Providing file pointers
^^^^^^^^^^^^^^^^^^^^^^^

There are a few functions in the SUNDIALS C API which take a ``FILE*`` argument.
Since there is no portable way to convert between a Fortran file descriptor and
a C file pointer, SUNDIALS provides two utility functions for creating a
``FILE*`` and destroying it. These functions are defined in the module
``fsundials_futils_mod``.

.. c:function:: FILE* SUNDIALSFileOpen(filename, mode)

   The function allocates a ``FILE*`` by calling the C function ``fopen`` with
   the provided filename and I/O mode.

   **Arguments:**
      * ``filename`` -- the full path to the file, that should have Fortran
        type ``character(kind=C_CHAR, len=*)``.
      * ``mode`` -- the I/O mode to use for the file.  This should have the
        Fortran type ``character(kind=C_CHAR, len=*)``.  The string begins
        with one of the following characters:

        * ``r``  to open a text file for reading
        * ``r+`` to open a text file for reading/writing
        * ``w``  to truncate a text file to zero length or create it for writing
        * ``w+`` to open a text file for reading/writing or create it if it does
          not exist
        * ``a``  to open a text file for appending, see documentation of ``fopen`` for
          your system/compiler
        * ``a+`` to open a text file for reading/appending, see documentation for
          ``fopen`` for your system/compiler

   **Return value:**
      * The function returns a ``type(C_PTR)`` which holds a C ``FILE*``.


.. c:function:: void SUNDIALSFileClose(fp)

   The function deallocates a C ``FILE*`` by calling the C function ``fclose``
   with the provided pointer.

   **Arguments:**
      * ``fp`` -- the C ``FILE*`` that was previously obtained from ``fopen``.
        This should have the Fortran type ``type(c_ptr)``.


.. _SUNDIALS.Fortran.Portability:

Important notes on portability
------------------------------

The SUNDIALS Fortran 2003 interface *should* be compatible with any compiler
supporting the Fortran 2003 ISO standard. However, it has only been tested and
confirmed to be working with GNU Fortran 4.9+ and Intel Fortran 18.0.1+.

Upon compilation of SUNDIALS, Fortran module (``.mod``) files are generated for
each Fortran 2003 interface. These files are highly compiler specific, and thus
it is almost always necessary to compile a consuming application with the same
compiler that was used to generate the modules.


.. _SUNDIALS.Fortran.CommonIssues:

Common Issues
-------------

In this subsection, we list some common issues users run into when using the Fortran
interfaces.


**Strange Segmentation Fault in User-Supplied Functions**

One common issue we have seen trip up users (and even ourselves) has the symptom
of segmentation fault in a user-supplied function (such as the RHS) when trying
to use one of the callback arguments. For example, in the following RHS
function, we will get a segfault on line 21:

.. code-block:: fortran
   :linenos:
   :emphasize-lines: 8, 21

   integer(c_int) function ff(t, yvec, ydotvec, user_data) &
      result(ierr) bind(C)

      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      implicit none

      real(c_double) :: t ! <===== Missing value attribute
      type(N_Vector) :: yvec
      type(N_Vector) :: ydotvec
      type(c_ptr)    :: user_data

      real(c_double) :: e
      real(c_double) :: u, v
      real(c_double) :: tmp1, tmp2
      real(c_double), pointer :: yarr(:)
      real(c_double), pointer :: ydotarr(:)

      ! get N_Vector data arrays
      yarr => FN_VGetArrayPointer(yvec)
      ydotarr => FN_VGetArrayPointer(ydotvec) ! <===== SEGFAULTS HERE

      ! extract variables
      u = yarr(1)
      v = yarr(2)

      ! fill in the RHS function:
      !  [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
      !  [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))]
      tmp1 = (-ONE+u*u-r(t))/(TWO*u)
      tmp2 = (-TWO+v*v-s(t))/(TWO*v)
      ydotarr(1) = ZERO
      ydotarr(2) = e*tmp1 - tmp2 + sdot(t)/(TWO*vtrue(t))

      ! return success
      ierr = 0
      return

   end function


The subtle bug in the code causing the segfault is on line 8. It should read
``real(c_double), value :: t`` instead of ``real(c_double) :: t`` (notice the
``value`` attribute). Fundamental types that are passed by value in C need
the ``value`` attribute.
