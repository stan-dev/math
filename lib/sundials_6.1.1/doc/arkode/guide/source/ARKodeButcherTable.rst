.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKodeButcherTable:

==============================
Butcher Table Data Structure
==============================

To store the Butcher table defining a Runge--Kutta method ARKODE provides the
:c:type:`ARKodeButcherTable` type and several related utility routines. We use
the following Butcher table notation (shown for a 3-stage method):

.. math::

   \begin{array}{r|c}
     c & A \\
     \hline
     q & b \\
     p & \tilde{b}
   \end{array}
   \quad = \quad
   \begin{array}{r|ccc}
     c_1 & a_{1,1} & a_{1,2} & a_{1,3} \\
     c_2 & a_{2,1} & a_{2,2} & a_{2,3} \\
     c_3 & a_{3,1} & a_{3,2} & a_{3,3} \\
     \hline
     q & b_1 & b_2 & b_3 \\
     p & \tilde{b}_1 & \tilde{b}_2 & \tilde{b}_3
   \end{array}

where the method and embedding share stage :math:`A` and abscissa :math:`c`
values, but use their stages :math:`z_i` differently through the coefficients
:math:`b` and :math:`\tilde{b}` to generate methods of orders :math:`q` (the
main method) and :math:`p` (the embedding, typically :math:`q = p+1`, though
sometimes this is reversed). :c:type:`ARKodeButcherTable` is defined as

.. c:type:: ARKodeButcherTableMem* ARKodeButcherTable

where ``ARKodeButcherTableMem`` is the structure

.. code-block:: c

   typedef struct ARKodeButcherTableMem {

     int q;
     int p;
     int stages;
     realtype **A;
     realtype *c;
     realtype *b;
     realtype *d;

   };

where ``stages`` is the number of stages in the RK method, the variables ``q``,
``p``, ``A``, ``c``, and ``b`` have the same meaning as in the Butcher table
above, and ``d`` is used to store :math:`\tilde{b}`.

.. _ARKodeButcherTable.Functions:

ARKodeButcherTable functions
-----------------------------

.. _ARKodeButcherTable.FunctionsTable:
.. table:: ARKodeButcherTable functions

   +----------------------------------------------+------------------------------------------------------------+
   | **Function name**                            | **Description**                                            |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_LoadERK()`       | Retrieve a given explicit Butcher table by its unique name |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_LoadDIRK()`      | Retrieve a given implicit Butcher table by its unique name |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Alloc()`         | Allocate an empty Butcher table                            |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Create()`        | Create a new Butcher table                                 |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Copy()`          | Create a copy of a Butcher table                           |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Space()`         | Get the Butcher table real and integer workspace size      |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Free()`          | Deallocate a Butcher table                                 |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_Write()`         | Write the Butcher table to an output file                  |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_CheckOrder()`    | Check the order of a Butcher table                         |
   +----------------------------------------------+------------------------------------------------------------+
   | :c:func:`ARKodeButcherTable_CheckARKOrder()` | Check the order of an ARK pair of Butcher tables           |
   +----------------------------------------------+------------------------------------------------------------+

.. c:function:: ARKodeButcherTable ARKodeButcherTable_LoadERK(ARKODE_ERKTableID emethod)

   Retrieves a specified explicit Butcher table. The prototype for this
   function, as well as the integer names for each provided method, are defined
   in the header file ``arkode/arkode_butcher_erk.h``.  For further information
   on these tables and their corresponding identifiers, see :numref:`Butcher`.

   **Arguments:**
      * *emethod* -- integer input specifying the given Butcher table.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *emethod* was invalid.


.. c:function:: ARKodeButcherTable ARKodeButcherTable_LoadDIRK(ARKODE_DIRKTableID imethod)

   Retrieves a specified diagonally-implicit Butcher table. The prototype for
   this function, as well as the integer names for each provided method, are
   defined in the header file ``arkode/arkode_butcher_dirk.h``.  For further
   information on these tables and their corresponding identifiers, see
   :numref:`Butcher`.

   **Arguments:**
      * *imethod* -- integer input specifying the given Butcher table.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *imethod* was invalid.

.. c:function:: ARKodeButcherTable ARKodeButcherTable_Alloc(int stages, booleantype embedded)

   Allocates an empty Butcher table.

   **Arguments:**
      * *stages* -- the number of stages in the Butcher table.
      * *embedded* -- flag denoting whether the Butcher table has an embedding
        (``SUNTRUE``) or not (``SUNFALSE``).

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *stages* was invalid or an allocation error occurred.

.. c:function:: ARKodeButcherTable ARKodeButcherTable_Create(int s, int q, int p, realtype *c, realtype *A, realtype *b, realtype *d)

   Allocates a Butcher table and fills it with the given values.

   **Arguments:**
      * *s* -- number of stages in the RK method.
      * *q* -- global order of accuracy for the RK method.
      * *p* -- global order of accuracy for the embedded RK method.
      * *c* -- array (of length *s*) of stage times for the RK method.
      * *A* -- array of coefficients defining the RK stages. This should be
        stored as a 1D array of size *s*s*, in row-major order.
      * *b* -- array of coefficients (of length *s*) defining the time step solution.
      * *d* -- array of coefficients (of length *s*) defining the embedded solution.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer if *stages* was invalid or an allocation error occurred.

   **Notes:**
      If the method does not have an embedding then *d* should be
      ``NULL`` and *p* should be equal to zero.

      .. warning::
         When calling this function from Fortran, it is important to note that ``A`` is expected
         to be in row-major ordering.

.. c:function:: ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B)

   Creates copy of the given Butcher table.

   **Arguments:**
      * *B* -- the Butcher table to copy.

   **Return value:**
      * :c:type:`ARKodeButcherTable` structure if successful.
      * ``NULL`` pointer an allocation error occurred.

.. c:function:: void ARKodeButcherTable_Space(ARKodeButcherTable B, sunindextype *liw, sunindextype *lrw)

   Get the real and integer workspace size for a Butcher table.

   **Arguments:**
      * *B* -- the Butcher table.
      * *lenrw* -- the number of ``realtype`` values in the Butcher table workspace.
      * *leniw* -- the number of integer values in the Butcher table workspace.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_MEM_NULL* if the Butcher table memory was ``NULL``.

.. c:function:: void ARKodeButcherTable_Free(ARKodeButcherTable B)

   Deallocate the Butcher table memory.

   **Arguments:**
      * *B* -- the Butcher table.

.. c:function:: void ARKodeButcherTable_Write(ARKodeButcherTable B, FILE *outfile)

   Write the Butcher table to the provided file pointer.

   **Arguments:**
      * *B* -- the Butcher table.
      * *outfile* -- pointer to use for printing the Butcher table.

   **Notes:**
      The *outfile* argument can be ``stdout`` or ``stderr``, or it
      may point to a specific file created using ``fopen``.

.. c:function:: int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int* q, int* p, FILE* outfile)

   Determine the analytic order of accuracy for the specified Butcher
   table. The analytic (necessary) conditions are checked up to order 6. For
   orders greater than 6 the Butcher simplifying (sufficient) assumptions are
   used.

   **Arguments:**
      * *B* -- the Butcher table.
      * *q* -- the measured order of accuracy for the method.
      * *p* -- the measured order of accuracy for the embedding; 0 if the
        method does not have an embedding.
      * *outfile* -- file pointer for printing results; ``NULL`` to suppress
        output.

   **Return value:**
      * *0* -- success, the measured vales of *q* and *p* match the values of
        *q* and *p* in the provided Butcher tables.
      * *1* -- warning, the values of *q* and *p* in the provided Butcher tables
        are *lower* than the measured values, or the measured values achieve the
        *maximum order* possible with this function and the values of *q* and
        *p* in the provided Butcher tables table are higher.
      * *-1* -- failure, the values of *q* and *p* in the provided Butcher tables
        are *higher* than the measured values.
      * *-2* -- failure, the input Butcher table or critical table contents are
        ``NULL``.

   **Notes:**
      For embedded methods, if the return flags for *q* and *p* would
      differ, failure takes precedence over warning, which takes precedence over
      success.


.. c:function:: int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1, ARKodeButcherTable B2, int *q, int *p, FILE *outfile)

   Determine the analytic order of accuracy (up to order 6) for a specified
   ARK pair of Butcher tables.

   **Arguments:**
      * *B1* -- a Butcher table in the ARK pair.
      * *B2* -- a Butcher table in the ARK pair.
      * *q* -- the measured order of accuracy for the method.
      * *p* -- the measured order of accuracy for the embedding; 0 if the
        method does not have an embedding.
      * *outfile* -- file pointer for printing results; ``NULL`` to suppress
        output.

   **Return value:**
      * *0* -- success, the measured vales of *q* and *p* match the values of
        *q* and *p* in the provided Butcher tables.
      * *1* -- warning, the values of *q* and *p* in the provided Butcher tables
        are *lower* than the measured values, or the measured values achieve the
        *maximum order* possible with this function and the values of *q* and
        *p* in the provided Butcher tables table are higher.
      * *-1* -- failure, the input Butcher tables or critical table contents are
        ``NULL``.

   **Notes:**
      For embedded methods, if the return flags for *q* and *p* would
      differ, warning takes precedence over success.
