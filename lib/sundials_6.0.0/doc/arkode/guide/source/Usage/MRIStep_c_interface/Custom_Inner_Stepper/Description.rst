.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Description:

The MRIStepInnerStepper Class
-----------------------------

As with other SUNDIALS classes, the :c:type:`MRIStepInnerStepper` abstract base
class is implemented using a C structure containing a ``content`` pointer to the
derived class member data and a structure of function pointers the derived class
implementations of the virtual methods. The :c:type:`MRIStepInnerStepper`
type is defined in ``include/arkode/arkode.h`` as

.. c:type:: struct _MRIStepInnerStepper *MRIStepInnerStepper

The actual definitions of the ``_MRIStepInnerStepper`` structure and the
corresponding operations structure are kept private to allow for the object
internals to change without impacting user code. The following sections describe
the :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods` and the virtual
:numref:`ARKODE.Usage.MRIStep.CustomInnerStepper.Description.ImplMethods` that a must be
provided by a derived class.

.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods:

Base Class Methods
^^^^^^^^^^^^^^^^^^

This section describes methods provided by the :c:type:`MRIStepInnerStepper`
abstract base class that aid the user in implementing derived classes. This
includes functions for creating and destroying a generic base class object,
attaching and retrieving the derived class ``content`` pointer, setting function
pointers to derived class method implementations, and accessing base class data
e.g., for computing the forcing term :eq:`ARKODE_MRI_forcing_poly`.

.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods.CreateDestroy:

Creating and Destroying an Object
"""""""""""""""""""""""""""""""""

.. c:function:: int MRIStepInnerStepper_Create(SUNContext sunctx, MRIStepInnerStepper *stepper)

   This function creates an :c:type:`MRIStepInnerStepper` object to which a user
   should attach the member data (content) pointer and method function pointers.

   **Arguments:**
      * ``sunctx`` -- the SUNDIALS simulation context.
      * ``stepper`` -- a pointer to an inner stepper object.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_MEM_FAIL if a memory allocation error occurs

   **Example usage:**

   .. code-block:: C

      /* create an instance of the base class */
      MRIStepInnerStepper inner_stepper = NULL;
      flag = MRIStepInnerStepper_Create(&inner_stepper);

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``

   .. note::

      See :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods.Content` and
      :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods.AttachFunctions`
      for details on how to attach member data and method function pointers.


.. c:function:: int MRIStepInnerStepper_Free(MRIStepInnerStepper *stepper)

   This function destroys an :c:type:`MRIStepInnerStepper` object.

   **Arguments:**
      * *stepper* -- a pointer to an inner stepper object.

   **Return value:**
      * ARK_SUCCESS if successful

   **Example usage:**

   .. code-block:: C

      /* destroy an instance of the base class */
      flag = MRIStepInnerStepper_Free(&inner_stepper);

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``

   .. note::

      This function only frees memory allocated within the base class and the
      base class structure itself. The user is responsible for freeing any
      memory allocated for the member data (content).

.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods.Content:

Attaching and Accessing the Content Pointer
"""""""""""""""""""""""""""""""""""""""""""

.. c:function:: int MRIStepInnerStepper_SetContent(MRIStepInnerStepper stepper, void *content)

   This function attaches a member data (content) pointer to an
   :c:type:`MRIStepInnerStepper` object.

   **Arguments:**
      * *stepper* -- an inner stepper object.
      * *content* -- a pointer to the stepper member data.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_ILL_INPUT if the stepper is ``NULL``

   **Example usage:**

   .. code-block:: C

      /* set the inner stepper content pointer */
      MyStepperContent my_object_data;
      flag = MRIStepInnerStepper_SetContent(inner_stepper, &my_object_data);

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. c:function:: int MRIStepInnerStepper_GetContent(MRIStepInnerStepper stepper, void **content)

   This function retrieves the member data (content) pointer from an
   :c:type:`MRIStepInnerStepper` object.

   **Arguments:**
      * *stepper* -- an inner stepper object.
      * *content* -- a pointer to set to the stepper member data pointer.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_ILL_INPUT if the stepper is ``NULL``

   **Example usage:**

   .. code-block:: C

      /* get the inner stepper content pointer */
      void             *content;
      MyStepperContent *my_object_data;

      flag = MRIStepInnerStepper_GetContent(inner_stepper, &content);
      my_object_data = (MyStepperContent*) content;

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods.AttachFunctions:

Setting Member Functions
""""""""""""""""""""""""

.. c:function:: int MRIStepInnerStepper_SetEvolveFn(MRIStepInnerStepper stepper, MRIStepInnerEvolveFn fn)

   This function attaches an :c:type:`MRIStepInnerEvolveFn` function to an
   :c:type:`MRIStepInnerStepper` object.

   **Arguments:**
      * *stepper* -- an inner stepper object.
      * *fn* -- the :c:type:`MRIStepInnerStepper` function to attach.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_ILL_INPUT if the stepper is ``NULL``

   **Example usage:**

   .. code-block:: C

      /* set the inner stepper evolve function */
      flag = MRIStepInnerStepper_SetEvolveFn(inner_stepper, MyEvolve);

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. c:function:: int MRIStepInnerStepper_SetFullRhsFn(MRIStepInnerStepper stepper, MRIStepInnerFullRhsFn fn)

   This function attaches an :c:type:`MRIStepInnerFullRhsFn` function to an
   :c:type:`MRIStepInnerStepper` object.

   **Arguments:**
      * *stepper* -- an inner stepper object.
      * *fn* -- the :c:type:`MRIStepInnerFullRhsFn` function to attach.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_ILL_INPUT if the stepper is ``NULL``

   **Example usage:**

   .. code-block:: C

      /* set the inner stepper full right-hand side function */
      flag = MRIStepInnerStepper_SetFullRhsFn(inner_stepper, MyFullRHS);

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. c:function:: int MRIStepInnerStepper_SetResetFn(MRIStepInnerStepper stepper, MRIStepInnerResetFn fn)

   This function attaches an :c:type:`MRIStepInnerResetFn` function to an
   :c:type:`MRIStepInnerStepper` object.

   **Arguments:**
      * *stepper* -- an inner stepper object.
      * *fn* -- the :c:type:`MRIStepInnerResetFn` function to attach.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_ILL_INPUT if the stepper is ``NULL``

   **Example usage:**

   .. code-block:: C

      /* set the inner stepper reset function */
      flag = MRIStepInnerStepper_SetResetFn(inner_stepper, MyReset);

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``

.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods.Forcing:

Applying and Accessing Forcing Data
"""""""""""""""""""""""""""""""""""

When integrating the ODE :eq:`ARKODE_MRI_IVP` the :c:type:`MRIStepInnerStepper` is
responsible for evaluating ODE right-hand side function :math:`f^F(t,v)` as well
as computing and applying the forcing term :eq:`ARKODE_MRI_forcing_poly` to obtain the
full right-hand side of the inner (fast) ODE :eq:`ARKODE_MRI_IVP`. The functions in
this section can be used to either apply the inner (fast) forcing or access the
data necessary to construct the inner (fast) forcing polynomial.


.. c:function:: int MRIStepInnerStepper_AddForcing(MRIStepInnerStepper stepper, realtype t, N_Vector ff)

   This function computes the forcing term :eq:`ARKODE_MRI_forcing_poly` at the input
   time *t* and adds it to input vector *ff*, i.e., the inner (fast) right-hand
   side vector.

   **Arguments:**
      * *stepper* -- an inner stepper object.
      * *t* -- the time at which the forcing should be evaluated.
      * *f* -- the vector to which the forcing should be applied.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_ILL_INPUT if the stepper is ``NULL``

   **Example usage:**

   .. code-block:: C

      /* compute the forcing term and add it the fast RHS vector */
      flag = MRIStepInnerStepper_AddForcing(inner_stepper, t, f_fast);

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. c:function:: int MRIStepInnerStepper_GetForcingData(MRIStepInnerStepper stepper, realtype *tshift, realtype *tscale, N_Vector **forcing, int *nforcing)

   This function provides access to data necessary to compute the forcing term
   :eq:`ARKODE_MRI_forcing_poly`. This includes the shift and scaling factors for the
   normalized time :math:`\tau = (t - t_{n,i-1}^S)/(h^S \Delta c_i^S)` and the
   array of polynomial coefficient vectors :math:`\hat{\gamma}^{\{k\}}_i`.

   **Arguments:**
      * *stepper* -- an inner stepper object.
      * *tshift* -- the time shift to apply to the current time when computing the
        forcing, :math:`t_{n,i-1}^S`.
      * *tscale* -- the time scaling to apply to the current time when computing
        the forcing, :math:`h^S \Delta c_i^S`.
      * *forcing* -- a pointer to an array of forcing vectors,
        :math:`\hat{\gamma}^{\{k\}}_i`.
      * *nforcing* -- the number of forcing vectors.

   **Return value:**
      * ARK_SUCCESS if successful
      * ARK_ILL_INPUT if the stepper is ``NULL``

   **Example usage:**

   .. code-block:: C

      int      k, flag;
      int      nforcing_vecs;   /* number of forcing vectors */
      double   tshift, tscale;  /* time normalization values */
      double   tau;             /* normalized time           */
      double   tau_k;           /* tau raised to the power k */
      N_Vector *forcing_vecs;   /* array of forcing vectors  */

      /* get the forcing data from the inner (fast) stepper */
      flag = MRIStepInnerStepper_GetForcingData(inner_stepper, &tshift, &tscale,
                                                &forcing_vecs, &nforcing_vecs);

      /* compute the normalized time, initialize tau^k */
      tau   = (t - tshift) / tscale;
      tau_k = 1.0;

      /* compute the polynomial forcing terms and add them to fast RHS vector */
      for (k = 0; k < nforcing_vecs; k++)
      {
        N_VLinearSum(1.0, f_fast, tau_k, forcing_vecs[k], f_fast);
        tau_k *= tau;
      }

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Description.ImplMethods:

Implementation Specific Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the required and optional virtual methods defined by the
:c:type:`MRIStepInnerStepper` abstract base class.

Required Member Functions
"""""""""""""""""""""""""

An :c:type:`MRIStepInnerStepper` *must* provide implementations of the following
member functions:


.. c:type:: int (*MRIStepInnerEvolveFn)(MRIStepInnerStepper stepper, realtype t0, realtype tout, N_Vector v)

   This function advances the state vector *v* for the inner (fast) ODE system
   from time *t0* to time *tout*.

   **Arguments:**
      * *stepper* -- the inner stepper object.
      * *t0* -- the initial time for the inner (fast) integration.
      * *tout* -- the final time for the inner (fast) integration.
      * *v* -- on input the state at time *t0* and, on output, the state at time
        *tout*.

   **Return value:**
      An :c:type:`MRIStepInnerEvolveFn` should return 0 if successful, a positive
      value if a recoverable error occurred, or a negative value if it failed
      unrecoverably.

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. c:type:: int (*MRIStepInnerFullRhsFn)(MRIStepInnerStepper stepper, realtype t, N_Vector v, N_Vector f, int mode)

   This function computes the full right-hand side function of the inner (fast)
   ODE, :math:`f^F(t,v)` in :eq:`ARKODE_MRI_IVP` for a given value of the independent
   variable *t* and state vector *y*.

   **Arguments:**
      * *stepper* -- the inner stepper object.
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *f* -- the output vector that forms a portion the ODE right-hand side,
        :math:`f^F(t,y)` in :eq:`ARKODE_IVP_two_rate`.
      * *mode* -- a flag indicating the purpose for which the right-hand side
        function evaluation is called.

        * ``ARK_FULLRHS_START`` -- called at the beginning of the simulation
        * ``ARK_FULLRHS_END``   -- called at the end of a successful step
        * ``ARK_FULLRHS_OTHER`` -- called elsewhere e.g., for dense output

   **Return value:**
      An :c:type:`MRIStepInnerFullRhsFn` should return 0 if successful, a positive
      value if a recoverable error occurred, or a negative value if it failed
      unrecoverably.

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``

Optional Member Functions
"""""""""""""""""""""""""

An :c:type:`MRIStepInnerStepper` *may* provide implementations of any of the
following member functions:

.. c:type:: int (*MRIStepInnerResetFn)(MRIStepInnerStepper stepper, realtype tR, N_Vector vR)

   This function resets the inner (fast) stepper state to the provided
   independent variable value and dependent variable vector.

   **Arguments:**
      * *stepper* -- the inner stepper object.
      * *tR* -- the value of the independent variable :math:`t_R`.
      * *vR* -- the value of the dependent variable vector :math:`v(t_R)`.

   **Return value:**
      An :c:type:`MRIStepInnerResetFn` should return 0 if successful, a positive
      value if a recoverable error occurred, or a negative value if it failed
      unrecoverably.

   **Example codes:**
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``
