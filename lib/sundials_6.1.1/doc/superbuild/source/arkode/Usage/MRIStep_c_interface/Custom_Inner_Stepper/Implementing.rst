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

.. _ARKODE.Usage.MRIStep.CustomInnerStepper.Implementing:

Implementing an MRIStepInnerStepper
-----------------------------------

To create an MRIStepInnerStepper implementation:

#. Define the stepper-specific content.

   This is typically a user-defined structure in C codes, a user-defined class
   or structure in C++ codes, or a user-defined module in Fortran codes. This
   content should hold any data necessary to perform the operations defined by
   the :c:type:`MRIStepInnerStepper` member functions.

#. Define implementations of the required member functions (see
   :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper.Description.ImplMethods`).

   These are typically user-defined functions in C, member functions of the
   user-defined structure or class in C++, or functions contained in the
   user-defined module in Fortran.

   Note that all member functions are passed the :c:type:`MRIStepInnerStepper`
   object and the stepper-specific content can, if necessary, be retrieved using
   :c:func:`MRIStepInnerStepper_GetContent`.

#. In the user code, before creating the MRIStep memory structure with
   :c:func:`MRIStepCreate`, do the following:

   #. Create an :c:type:`MRIStepInnerStepper` object with
      :c:func:`MRIStepInnerStepper_Create`.

   #. Attach a pointer to the stepper content to the
      :c:type:`MRIStepInnerStepper` object with
      :c:func:`MRIStepInnerStepper_SetContent` if necessary, e.g., when the
      content is a C structure.

   #. Attach the member function implementations using the functions described
      in :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper.Description.BaseMethods.AttachFunctions`.

#. Attach the :c:type:`MRIStepInnerStepper` object to the MRIStep memory
   structure with :c:func:`MRIStepCreate`.

For an example of creating and attaching a user-defined inner stepper see
the example code ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``
where CVODE is wrapped as an :c:type:`MRIStepInnerStepper`.
