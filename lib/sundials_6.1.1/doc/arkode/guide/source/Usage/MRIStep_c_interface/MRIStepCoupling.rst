.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.MRIStep.MRIStepCoupling:

MRI Coupling Coefficients Data Structure
----------------------------------------

MRIStep supplies several built-in MIS, MRI-GARK, and IMEX-MRI-GARK methods, see
:numref:`ARKODE.Usage.MRIStep.MRIStepCoupling.Tables` for the current set of coupling
tables and their corresponding identifiers. Additionally, a user may supply a
custom set of slow-to-fast time scale coupling coefficients by constructing a
coupling table and attaching it with :c:func:`MRIStepSetCoupling()`.

As described in :numref:`ARKODE.Mathematics.MRIStep`, the coupling from the slow time
scale to the fast time scale is encoded by a vector of slow 'stage time'
abscissae, :math:`c^S \in \mathbb{R}^{s+1}` and a set of coupling matrices
:math:`\Gamma^{\{k\}}\in\mathbb{R}^{(s+1)\times(s+1)}` and
:math:`\Omega^{\{k\}}\in\mathbb{R}^{(s+1)\times(s+1)}`. An ``MRIStepCoupling``
object stores this information and provides several related utility functions
for creating a coupling table. The ``MRIStepCoupling`` type is defined as:

.. c:type:: MRIStepCouplingMem *MRIStepCoupling

where ``MRIStepCouplingMem`` is the structure

.. code-block:: c

   struct MRIStepCouplingMem
   {
     int nmat;
     int stages;
     int q;
     int p;
     realtype ***G;
     realtype ***W;
     realtype *c;
   };

and the members of the strucutre are:

   * ``nmat`` corresponds to the number of coupling matrices
     :math:`\Omega^{\{k\}}` for the slow-nonstiff terms and/or
     :math:`\Gamma^{\{k\}}` for the slow-stiff terms in :eq:`ARKODE_IVP_two_rate`,

   * ``stages`` is the number of abscissae i.e., :math:`s+1` above,

   * ``q`` and ``p`` indicate the orders of accuracy for both the method and
     the embedding, respectively,

   * ``W`` is a three-dimensional array with dimensions
     ``[nmat][stages][stages]`` containing the method's :math:`\Omega^{\{k\}}`
     coupling matrices for the slow-nonstiff (explicit) terms in
     :eq:`ARKODE_IVP_two_rate`,

   * ``G`` is a three-dimensional array with dimensions
     ``[nmat][stages][stages]`` containing the method's :math:`\Gamma^{\{k\}}`
     coupling matrices for the slow-stiff (implicit) terms in
     :eq:`ARKODE_IVP_two_rate`, and

   * ``c`` is an array of length ``stages`` containing the slow abscissae
     :math:`c^S` for the method.


.. _ARKODE.Usage.MRIStep.MRIStepCoupling.Functions:

MRIStepCoupling functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the functions for creating and interacting with coupling
tables. The function prototypes and as well as the relevant integer constants
are defined ``arkode/arkode_mristep.h``.

.. _ARKODE.Usage.MRIStep.MRIStepCoupling.Functions.Table:
.. table:: MRIStepCoupling functions

   +---------------------------------------+--------------------------------------------------------------------+
   | Function name                         | Description                                                        |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_LoadTable()` | Loads a pre-defined MRIStepCoupling table                          |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Alloc()`     | Allocate an empty MRIStepCoupling table                            |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Create()`    | Create a new MRIStepCoupling table from coefficients               |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_MIStoMRI()`  | Create a new MRIStepCoupling table from a slow Butcher table       |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Copy()`      | Create a copy of a MRIStepCoupling table                           |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Space()`     | Get the MRIStepCoupling table real and integer workspace sizes     |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Free()`      | Deallocate a MRIStepCoupling table                                 |
   +---------------------------------------+--------------------------------------------------------------------+
   | :c:func:`MRIStepCoupling_Write()`     | Write the MRIStepCoupling table to an output file                  |
   +---------------------------------------+--------------------------------------------------------------------+


.. c:function:: MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID imethod)

   Retrieves a specified coupling table. For further information on the current
   set of coupling tables and their corresponding identifiers, see
   :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling.Tables`.


   **Arguments:**
      * ``itable`` -- the coupling table identifier.

   **Return value:**
      * An :c:type:`MRIStepCoupling` structure if successful.
      * A ``NULL`` pointer if *itable* was invalid or an allocation error occurred.


.. c:function:: MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages, int type)

   Allocates an empty MRIStepCoupling table.

   **Arguments:**
      * ``nmat`` -- number of :math:`\Omega^{\{k\}}` and/or :math:`\Gamma^{\{k\}}`
        matrices in the coupling table.
      * ``stages`` -- number of stages in the coupling table.
      * ``type`` -- the method type: explicit (0), implicit (1), or ImEx (2).

   **Return value:**
      * An :c:type:`MRIStepCoupling` structure if successful.
      * A ``NULL`` pointer if *stages* or *type* was invalid or an allocation error
        occurred.

   .. note::

      For explicit methods only the W array is allocated, with implicit methods
      only the G array is allocated, and for ImEx methods both W and G are
      allocated.


.. c:function:: MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages, int q, int p, realtype *W, realtype *G, realtype *c)

   Allocates a coupling table and fills it with the given values.

   **Arguments:**
      * ``nmat`` -- number of :math:`\Omega^{\{k\}}` and/or :math:`\Gamma^{\{k\}}`
        matrices in the coupling table.
      * ``stages`` -- number of stages in the method.
      * ``q`` -- global order of accuracy for the method.
      * ``p`` -- global order of accuracy for the embedded method.
      * ``W`` -- array of coefficients defining the explicit coupling matrices
        :math:`\Omega^{\{k\}}`. The entries should be stored as a 1D array of size
        ``nmat * stages * stages``, in row-major order. If the slow method is
        implicit pass ``NULL``.
      * ``G`` -- array of coefficients defining the implicit coupling matrices
        :math:`\Gamma^{\{k\}}`. The entries should be stored as a 1D array of size
        ``nmat * stages * stages``, in row-major order. If the slow method is
        explicit pass ``NULL``.
      * ``c`` -- array of slow abscissae for the MRI method. The entries should be
        stored as a 1D array of length ``stages``.

   **Return value:**
      * An :c:type:`MRIStepCoupling` structure if successful.
      * A ``NULL`` pointer if ``stages`` was invalid, an allocation error occurred,
        or the input data arrays are inconsistent with the method type.

   .. note::

      As embeddings are not currently supported in MRIStep, ``p`` should be
      equal to zero.

.. c:function:: MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B, int q, int p)

   Creates an MRI coupling table for a traditional MIS method based on the slow
   Butcher table *B*, following the formula shown in :eq:`ARKODE_MIS_to_MRI`

   **Arguments:**
      * ``B`` -- the :c:type:`ARKodeButcherTable` for the 'slow' MIS method.
      * ``q`` -- the overall order of the MIS/MRI method.
      * ``p`` -- the overall order of the MIS/MRI embedding.

   **Return value:**
      * An :c:type:`MRIStepCoupling` structure if successful.
      * A ``NULL`` pointer if an allocation error occurred.

   .. note::

      The :math:`s`-stage slow Butcher table must have an explicit first stage
      (i.e., :math:`c_1=0` and :math:`A_{1,j}=0` for :math:`1\le j\le s`) and
      sorted abscissae (i.e., :math:`c_{i} \ge c_{i-1}` for :math:`2\le i\le s`).

      Since an MIS method is at most third order accurate, and even then only if
      it meets certain compatibility criteria (see :eq:`ARKODE_MIS_order3`), the values
      of *q* and *p* may differ from the method and embedding orders of accuracy
      for the Runge--Kutta method encoded in *B*, which is why these arguments
      should be supplied separately.

      As embeddings are not currently supported in MRIStep, then *p* should be
      equal to zero.


.. c:function:: MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling C)

   Creates copy of the given coupling table.

   **Arguments:**
      * ``C`` -- the coupling table to copy.

   **Return value:**
      * An :c:type:`MRIStepCoupling` structure if successful.
      * A ``NULL`` pointer if an allocation error occurred.


.. c:function:: void MRIStepCoupling_Space(MRIStepCoupling C, sunindextype *liw, sunindextype *lrw)

   Get the real and integer workspace size for a coupling table.

   **Arguments:**
      * ``C`` -- the coupling table.
      * ``lenrw`` -- the number of ``realtype`` values in the coupling table
        workspace.
      * ``leniw`` -- the number of integer values in the coupling table workspace.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_MEM_NULL* if the Butcher table memory was ``NULL``.


.. c:function:: void MRIStepCoupling_Free(MRIStepCoupling C)

   Deallocate the coupling table memory.

   **Arguments:**
      * ``C`` -- the coupling table.


.. c:function:: void MRIStepCoupling_Write(MRIStepCoupling C, FILE *outfile)

   Write the coupling table to the provided file pointer.

   **Arguments:**
      * ``C`` -- the coupling table.
      * ``outfile`` -- pointer to use for printing the table.

   .. note::

      The *outfile* argument can be ``stdout`` or ``stderr``, or it may point to
      a specific file created using ``fopen``.





.. _ARKODE.Usage.MRIStep.MRIStepCoupling.Tables:

MRI Coupling Tables
^^^^^^^^^^^^^^^^^^^

MRIStep currently includes three classes of coupling tables: those that encode
methods that are explicit at the slow time scale, those that are
diagonally-implicit and solve-decoupled at the slow time scale, and those that
encode methods with an implicit-explicit method at the slow time scale.  We list
the current identifiers, multirate order of accuracy, and relevant references
for each in the tables below. For methods with an implicit component, we also
list the number of implicit solves per step that are required at the slow time
scale.

Each of the coupling tables that are packaged with MRIStep are specified by a
unique ID having type:

.. c:type:: int ARKODE_MRITableID

with values specified for each method below (e.g., ``ARKODE_MIS_KW3``).



.. table:: Explicit MRI-GARK coupling tables. The default method for each order
           is marked with an asterisk (:math:`^*`).

   ==========================  ===========  =====================
   Table name                  Order        Reference
   ==========================  ===========  =====================
   ``ARKODE_MIS_KW3``          :math:`3^*`  :cite:p:`Schlegel:09`
   ``ARKODE_MRI_GARK_ERK33a``  3            :cite:p:`Sandu:19`
   ``ARKODE_MRI_GARK_ERK45a``  :math:`4^*`  :cite:p:`Sandu:19`
   ==========================  ===========  =====================


.. table:: Diagonally-implicit, solve-decoupled MRI-GARK coupling tables. The
           default method for each order is marked with an asterisk
           (:math:`^*`).

   =============================  ===========  ===============  ==================
   Table name                     Order        Implicit Solves  Reference
   =============================  ===========  ===============  ==================
   ``ARKODE_MRI_GARK_IRK21a``     :math:`2^*`  1                :cite:p:`Sandu:19`
   ``ARKODE_MRI_GARK_ESDIRK34a``  :math:`3^*`  3                :cite:p:`Sandu:19`
   ``ARKODE_MRI_GARK_ESDIRK46a``  :math:`4^*`  5                :cite:p:`Sandu:19`
   =============================  ===========  ===============  ==================


.. table:: Diagonally-implicit, solve-decoupled IMEX-MRI-GARK coupling tables.
           The default method for each order is marked with an asterisk
           (:math:`^*`).

   ===========================  ===========  ===============  ===================
   Table name                   Order        Implicit Solves  Reference
   ===========================  ===========  ===============  ===================
   ``ARKODE_IMEX_MRI_GARK3a``   :math:`3^*`  2                :cite:p:`ChiRen:21`
   ``ARKODE_IMEX_MRI_GARK3b``   3            2                :cite:p:`ChiRen:21`
   ``ARKODE_IMEX_MRI_GARK4``    :math:`4^*`  5                :cite:p:`ChiRen:21`
   ===========================  ===========  ===============  ===================
