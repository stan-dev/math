.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS.GPU:


Features for GPU Accelerated Computing
=============================================

.. ifconfig:: package_name == 'kinsol'

   In this section, we introduce the SUNDIALS GPU programming model and
   highlight SUNDIALS GPU features. The model leverages the fact that all of the
   SUNDIALS packages interact with simulation data either through the shared
   vector, matrix, and solver APIs (see Chapters :numref:`NVectors`,
   :numref:`SUNMatrix`, and :numref:`SUNLinSol`) or through user-supplied
   callback functions. Thus, under the model, the overall structure of the
   user’s calling program, and the way users interact with the SUNDIALS packages
   is similar to using SUNDIALS in CPU-only environments.

.. ifconfig:: package_name != 'kinsol'

   In this section, we introduce the SUNDIALS GPU programming model and
   highlight SUNDIALS GPU features. The model leverages the fact that all of the
   SUNDIALS packages interact with simulation data either through the shared
   vector, matrix, and solver APIs (see Chapters :numref:`NVectors`,
   :numref:`SUNMatrix`, :numref:`SUNLinSol`, and :numref:`SUNNonlinSol`) or
   through user-supplied callback functions. Thus, under the model, the overall
   structure of the user’s calling program, and the way users interact with the
   SUNDIALS packages is similar to using SUNDIALS in CPU-only environments.


.. _SUNDIALS.GPU.Model:

SUNDIALS GPU Programming Model
------------------------------

As described in :cite:p:`balos2021enabling`, within the SUNDIALS GPU programming
model, all control logic executes on the CPU, and all simulation data resides
wherever the vector or matrix object dictates as long as SUNDIALS is in control
of the program. That is, SUNDIALS will not migrate data (explicitly) from one
memory space to another. Except in the most advanced use cases, it is safe to
assume that data is kept resident in the GPU-device memory space. The
consequence of this is that, when control is passed from the user’s calling
program to SUNDIALS, simulation data in vector or matrix objects must be
up-to-date in the device memory space. Similarly, when control is passed from
SUNDIALS to the user’s calling program, the user should assume that any
simulation data in vector and matrix objects are up-to-date in the device memory
space. To put it succinctly, *it is the responsibility of the user’s calling
program to manage data coherency between the CPU and GPU-device memory spaces*
unless unified virtual memory (UVM), also known as managed memory, is being
utilized.  Typically, the GPU-enabled SUNDIALS modules provide functions to copy
data from the host to the device and vice-versa as well as support for unmanaged
memory or UVM. In practical terms, the way SUNDIALS handles distinct host and
device memory spaces means that *users need to ensure that the user-supplied
functions, e.g. the right-hand side function, only operate on simulation data in
the device memory space* otherwise extra memory transfers will be required and
performance will suffer.  The exception to this rule is if some form of hybrid
data partitioning (achievable with the NVECTOR_MANYVECTOR, see
:numref:`NVectors.ManyVector`) is utilized.

.. ifconfig:: package_name == 'kinsol'

   SUNDIALS provides many native shared features and modules that are
   GPU-enabled.  Currently, these include the NVIDIA CUDA platform
   :cite:p:`cuda_site`, AMD ROCm/HIP :cite:p:`rocm_site`, and Intel oneAPI
   :cite:p:`oneAPI_site`. :numref:`Usage.GPU.NVectorTable`--:numref:`Usage.GPU.SUNLinSolTable`
   summarize the shared SUNDIALS modules that are GPU-enabled, what GPU
   programming environments they support, and what class of memory they support
   (unmanaged or UVM).  Users may also supply their own GPU-enabled
   :c:type:`N_Vector`, :c:type:`SUNMatrix`, or :c:type:`SUNLinearSolver`
   implementation, and the capabilties will be leveraged since SUNDIALS operates
   on data through these APIs.

.. ifconfig:: package_name != 'kinsol'

   SUNDIALS provides many native shared features and modules that are
   GPU-enabled.  Currently, these include the NVIDIA CUDA platform
   :cite:p:`cuda_site`, AMD ROCm/HIP :cite:p:`rocm_site`, and Intel oneAPI
   :cite:p:`oneAPI_site`. :numref:`Usage.GPU.NVectorTable`--:numref:`Usage.GPU.SUNNonlinSolTable`
   summarize the shared SUNDIALS modules that are GPU-enabled, what GPU
   programming environments they support, and what class of memory they support
   (unmanaged or UVM).  Users may also supply their own GPU-enabled
   :c:type:`N_Vector`, :c:type:`SUNMatrix`, :c:type:`SUNLinearSolver`, or
   :c:type:`SUNNonlinearSolver` implementation, and the capabilties will be
   leveraged since SUNDIALS operates on data through these APIs.

In addition, SUNDIALS provides a memory management helper module
(see :numref:`SUNMemory`) to support applications which implement their own
memory management or memory pooling.

.. _Usage.GPU.NVectorTable:
.. table:: List of SUNDIALS GPU-enabled :c:type:`N_Vector` Modules

   ==========================================================  ===========  ===========  ===========  ================  ===========
   Module                                                      CUDA         ROCm/HIP     oneAPI       Unmanaged Memory  UVM
   ==========================================================  ===========  ===========  ===========  ================  ===========
   :ref:`NVECTOR_CUDA <NVectors.CUDA>`                         X                                      X                 X
   :ref:`NVECTOR_HIP  <NVectors.HIP>`                          X            X                         X                 X
   :ref:`NVECTOR_RAJA <NVectors.RAJA>`                         X            X            X            X                 X
   :ref:`NVECTOR_SYCL <NVectors.SYCL>`                         X\ :sup:`3`  X\ :sup:`3`  X            X                 X
   :ref:`NVECTOR_OPENMPDEV <NVectors.OPENMPDEV>`               X            X\ :sup:`2`  X\ :sup:`2`  X
   ==========================================================  ===========  ===========  ===========  ================  ===========

.. _Usage.GPU.SUNMatrixTable:
.. table:: List of SUNDIALS GPU-enabled :c:type:`SUNMatrix` Modules

   ==========================================================  ===========  ===========  ===========  ================  ===========
   Module                                                      CUDA         ROCm/HIP     oneAPI       Unmanaged Memory  UVM
   ==========================================================  ===========  ===========  ===========  ================  ===========
   :ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>`              X                                      X                 X
   :ref:`SUNMATRIX_MAGMADENSE <SUNMatrix.MagmaDense>`          X            X                         X                 X
   :ref:`SUNMATRIX_ONEMKLDENSE <SUNMatrix.OneMklDense>`        X\ :sup:`3`  X\ :sup:`3`  X            X                 X
   ==========================================================  ===========  ===========  ===========  ================  ===========

.. _Usage.GPU.SUNLinSolTable:
.. table:: List of SUNDIALS GPU-enabled :c:type:`SUNLinearSolver` Modules

   ==========================================================  ===========  ===========  ===========  ================  ===========
   Module                                                      CUDA         ROCm/HIP     oneAPI       Unmanaged Memory  UVM
   ==========================================================  ===========  ===========  ===========  ================  ===========
   :ref:`SUNLINSOL_CUSOLVERSP <SUNLinSol.cuSolverSp>`          X                                      X                 X
   :ref:`SUNLINSOL_MAGMADENSE <SUNLinSol.MagmaDense>`          X                                      X                 X
   :ref:`SUNLINSOL_ONEMKLDENSE <SUNLinSol.OneMklDense>`        X\ :sup:`3`  X\ :sup:`3`  X            X                 X
   :ref:`SUNLINSOL_SPGMR <SUNLinSol.SPGMR>`                    X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`       X\ :sup:`1`
   :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>`                  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`       X\ :sup:`1`
   :ref:`SUNLINSOL_SPTFQMR <SUNLinSol.SPTFQMR>`                X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`       X\ :sup:`1`
   :ref:`SUNLINSOL_SPBCGS <SUNLinSol.SPBCGS>`                  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`       X\ :sup:`1`
   :ref:`SUNLINSOL_PCG <SUNLinSol.PCG>`                        X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`       X\ :sup:`1`
   ==========================================================  ===========  ===========  ===========  ================  ===========

.. ifconfig:: package_name != 'kinsol'

   .. _Usage.GPU.SUNNonlinSolTable:
   .. table:: List of SUNDIALS GPU-enabled :c:type:`SUNNonlinearSolver` Modules

      ==========================================================  ===========  ===========  ===========  ================  ===========
      Module                                                      CUDA         ROCm/HIP     oneAPI       Unmanaged Memory  UVM
      ==========================================================  ===========  ===========  ===========  ================  ===========
      :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`            X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`       X\ :sup:`1`
      :ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>`    X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`  X\ :sup:`1`       X\ :sup:`1`
      ==========================================================  ===========  ===========  ===========  ================  ===========

Notes regarding the above tables:

1. This module inherits support from the NVECTOR module used
2. Support for ROCm/HIP and oneAPI are currently untested.
3. Support for CUDA and ROCm/HIP are currently untested.

In addition, note that implicit UVM (i.e. ``malloc`` returning UVM) is not
accounted for.


.. _SUNDIALS.GPU.Usage:

Steps for Using GPU Accelerated SUNDIALS
----------------------------------------

For any SUNDIALS package, the generalized steps a user needs to take to use GPU
accelerated SUNDIALS are:

#. Utilize a GPU-enabled ``N_Vector`` implementation. Initial data can be loaded
   on the host, but must be in the device memory space prior to handing control
   to SUNDIALS.

#. Utilize a GPU-enabled ``SUNLinearSolver`` linear solver (if applicable).

#. Utilize a GPU-enabled ``SUNMatrix`` implementation (if using a matrix-based
   linear solver).

#. Utilize a GPU-enabled ``SUNNonlinearSolver`` nonlinear solver (if
   applicable).

#. Write user-supplied functions so that they use data only in the device memory
   space (again, unless an atypical data partitioning is used). A few examples
   of these functions are the right-hand side evaluation function, the Jacobian
   evalution function, or the preconditioner evaulation function. In the context
   of CUDA and the right-hand side function, one way a user might ensure data is
   accessed on the device is, for example, calling a CUDA kernel, which does all
   of the computation, from a CPU function which simply extracts the underlying
   device data array from the :c:type:`N_Vector` object that is passed from
   SUNDIALS to the user-supplied function.

Users should refer to the above tables for a complete list of GPU-enabled
native SUNDIALS modules.
