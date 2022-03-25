..
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

.. _NVectors.Examples:

NVECTOR Examples
================

There are NVECTOR examples that may be installed for eac himplementation.
Each implementation makes use of the functions in ``test_nvector.c``.
These example functions show simple usage of the NVECTOR family
of functions. The input to the examples are the vector length, number
of threads (if threaded implementation), and a print timing flag.

The following is a list of the example functions in ``test_nvector.c``:

* ``Test_N_VClone``: Creates clone of vector and checks validity of clone.

* ``Test_N_VCloneEmpty``: Creates clone of empty vector and checks validity of clone.

* ``Test_N_VCloneVectorArray``: Creates clone of vector array and checks validity of cloned array.

* ``Test_N_VCloneVectorArray``: Creates clone of empty vector array and checks validity of cloned array.

* ``Test_N_VGetArrayPointer``: Get array pointer.

* ``Test_N_VSetArrayPointer``: Allocate new vector, set pointer to new vector array, and check values.

* ``Test_N_VGetLength``: Compares self-reported length to calculated length.

* ``Test_N_VGetCommunicator``: Compares self-reported communicator to the one used in constructor; or for MPI-unaware vectors it ensures that NULL is reported.

* ``Test_N_VLinearSum`` Case 1a: Test y =  x + y

* ``Test_N_VLinearSum`` Case 1b: Test y = -x + y

* ``Test_N_VLinearSum`` Case 1c: Test y = ax + y

* ``Test_N_VLinearSum`` Case 2a: Test x =  x + y

* ``Test_N_VLinearSum`` Case 2b: Test x =  x - y

* ``Test_N_VLinearSum`` Case 2c: Test x =  x + by

* ``Test_N_VLinearSum`` Case 3:  Test z =  x + y

* ``Test_N_VLinearSum`` Case 4a: Test z =  x - y

* ``Test_N_VLinearSum`` Case 4b: Test z = -x + y

* ``Test_N_VLinearSum`` Case 5a: Test z =  x + by

* ``Test_N_VLinearSum`` Case 5b: Test z = ax + y

* ``Test_N_VLinearSum`` Case 6a: Test z = -x + by

* ``Test_N_VLinearSum`` Case 6b: Test z = ax - y

* ``Test_N_VLinearSum`` Case 7:  Test z = a(x + y)

* ``Test_N_VLinearSum`` Case 8:  Test z = a(x - y)

* ``Test_N_VLinearSum`` Case 9:  Test z = ax + by

* ``Test_N_VConst``: Fill vector with constant and check result.

* ``Test_N_VProd``: Test vector multiply: z = x * y

* ``Test_N_VDiv``: Test vector division: z = x / y

* ``Test_N_VScale``: Case 1: scale: x = cx

* ``Test_N_VScale``: Case 2: copy: z = x

* ``Test_N_VScale``: Case 3: negate: z = -x

* ``Test_N_VScale``: Case 4: combination: z = cx

* ``Test_N_VAbs``: Create absolute value of vector.

* ``Test_N_VInv``: Compute z[i] = 1 / x[i]

** ``Test_N_VAddConst``: add constant vector: z = c + x

* ``Test_N_VDotProd``: Calculate dot product of two vectors.

* ``Test_N_VMaxNorm``: Create vector with known values, find and validate the max norm.

* ``Test_N_VWrmsNorm``: Create vector of known values, find and validate the weighted root mean square.

* ``Test_N_VWrmsNormMask``: Create vector of known values, find and validate the weighted root mean square using all elements except one.

* ``Test_N_VMin``: Create vector, find and validate the min.

* ``Test_N_VWL2Norm``: Create vector, find and validate the weighted Euclidean L2 norm.

* ``Test_N_VL1Norm``: Create vector, find and validate the L1 norm.

* ``Test_N_VCompare``: Compare vector with constant returning and validating comparison vector.

* ``Test_N_VInvTest``: Test z[i] = 1 / x[i]

* ``Test_N_VConstrMask``: Test mask of vector x with vector c.

* ``Test_N_VMinQuotient``: Fill two vectors with known values. Calculate and validate minimum quotient.

* ``Test_N_VLinearCombination``: Case 1a: Test x = a x

* ``Test_N_VLinearCombination``: Case 1b: Test z = a x

* ``Test_N_VLinearCombination``: Case 2a: Test x = a x + b y

* ``Test_N_VLinearCombination``: Case 2b: Test z = a x + b y

* ``Test_N_VLinearCombination``: Case 3a: Test x = x + a y + b z

* ``Test_N_VLinearCombination``: Case 3b: Test x = a x + b y + c z

* ``Test_N_VLinearCombination``: Case 3c: Test w = a x + b y + c z

* ``Test_N_VScaleAddMulti``: Case 1a: y = a x + y

* ``Test_N_VScaleAddMulti``: Case 1b: z = a x + y

* ``Test_N_VScaleAddMulti``: Case 2a: Y[i] = c[i] x + Y[i], i = 1,2,3

* ``Test_N_VScaleAddMulti``: Case 2b: Z[i] = c[i] x + Y[i], i = 1,2,3

* ``Test_N_VDotProdMulti``: Case 1: Calculate the dot product of two vectors

* ``Test_N_VDotProdMulti``: Case 2: Calculate the dot product of one vector with three other vectors in a vector array.

* ``Test_N_VLinearSumVectorArray``: Case 1: z = a x + b y

* ``Test_N_VLinearSumVectorArray``: Case 2a: Z[i] = a X[i] + b Y[i]

* ``Test_N_VLinearSumVectorArray``: Case 2b: X[i] = a X[i] + b Y[i]

* ``Test_N_VLinearSumVectorArray``: Case 2c: Y[i] = a X[i] + b Y[i]

* ``Test_N_VScaleVectorArray``: Case 1a: y = c y

* ``Test_N_VScaleVectorArray``: Case 1b: z = c y

* ``Test_N_VScaleVectorArray``: Case 2a: Y[i] = c[i] Y[i]

* ``Test_N_VScaleVectorArray``: Case 2b: Z[i] = c[i] Y[i]

* ``Test_N_VConstVectorArray``: Case 1a: z = c

* ``Test_N_VConstVectorArray``: Case 1b: Z[i] = c

* ``Test_N_VWrmsNormVectorArray``: Case 1a: Create a vector of know values, find and validate the weighted root mean square norm.

* ``Test_N_VWrmsNormVectorArray``: Case 1b: Create a vector array of three vectors of know values, find and validate the weighted root mean square norm of each.

* ``Test_N_VWrmsNormMaskVectorArray``: Case 1a: Create a vector of know values, find and validate the weighted root mean square norm using all elements except one.

* ``Test_N_VWrmsNormMaskVectorArray``: Case 1b: Create a vector array of three vectors of know values, find and validate the weighted root mean square norm of each using all elements except one.

* ``Test_N_VScaleAddMultiVectorArray``: Case 1a: y = a x + y

* ``Test_N_VScaleAddMultiVectorArray``: Case 1b: z = a x + y

* ``Test_N_VScaleAddMultiVectorArray``: Case 2a: Y[j][0] = a[j] X[0] + Y[j][0]

* ``Test_N_VScaleAddMultiVectorArray``: Case 2b: Z[j][0] = a[j] X[0] + Y[j][0]

* ``Test_N_VScaleAddMultiVectorArray``: Case 3a: Y[0][i] = a[0] X[i] + Y[0][i]

* ``Test_N_VScaleAddMultiVectorArray``: Case 3b: Z[0][i] = a[0] X[i] + Y[0][i]

* ``Test_N_VScaleAddMultiVectorArray``: Case 4a: Y[j][i] = a[j] X[i] + Y[j][i]

* ``Test_N_VScaleAddMultiVectorArray``: Case 4b: Z[j][i] = a[j] X[i] + Y[j][i]

* ``Test_N_VLinearCombinationVectorArray``: Case 1a: x = a x

* ``Test_N_VLinearCombinationVectorArray``: Case 1b: z = a x

* ``Test_N_VLinearCombinationVectorArray``: Case 2a: x = a x + b y

* ``Test_N_VLinearCombinationVectorArray``: Case 2b: z = a x + b y

* ``Test_N_VLinearCombinationVectorArray``: Case 3a: x = a x + b y + c z

* ``Test_N_VLinearCombinationVectorArray``: Case 3b: w = a x + b y + c z

* ``Test_N_VLinearCombinationVectorArray``: Case 4a: X[0][i] = c[0] X[0][i]

* ``Test_N_VLinearCombinationVectorArray``: Case 4b: Z[i] = c[0] X[0][i]

* ``Test_N_VLinearCombinationVectorArray``: Case 5a: X[0][i] = c[0] X[0][i] + c[1] X[1][i]

* ``Test_N_VLinearCombinationVectorArray``: Case 5b: Z[i] = c[0] X[0][i] + c[1] X[1][i]

* ``Test_N_VLinearCombinationVectorArray``: Case 6a: X[0][i] = X[0][i] + c[1] X[1][i] + c[2] X[2][i]

* ``Test_N_VLinearCombinationVectorArray``: Case 6b: X[0][i] = c[0] X[0][i] + c[1] X[1][i] + c[2] X[2][i]

* ``Test_N_VLinearCombinationVectorArray``: Case 6c: Z[i] = c[0] X[0][i] + c[1] X[1][i] + c[2] X[2][i]

* ``Test_N_VDotProdLocal``: Calculate MPI task-local portion of the dot product of two vectors.

* ``Test_N_VMaxNormLocal``: Create vector with known values, find and validate the MPI task-local portion of the max norm.

* ``Test_N_VMinLocal``: Create vector, find and validate the MPI task-local min.

* ``Test_N_VL1NormLocal``: Create vector, find and validate the MPI task-local portion of the L1 norm.

* ``Test_N_VWSqrSumLocal``: Create vector of known values, find and validate the MPI task-local portion of the weighted squared sum of two vectors.

* ``Test_N_VWSqrSumMaskLocal``: Create vector of known values, find and validate the MPI task-local portion of the weighted squared sum of two vectors, using all elements except one.

* ``Test_N_VInvTestLocal``: Test the MPI task-local portion of z[i] = 1 / x[i]

* ``Test_N_VConstrMaskLocal``: Test the MPI task-local portion of the mask of vector x with vector c.

* ``Test_N_VMinQuotientLocal``: Fill two vectors with known values. Calculate and validate the MPI task-local minimum quotient.

* ``Test_N_VMBufSize``: Tests for accuracy in the reported buffer size.

* ``Test_N_VMBufPack``: Tests for accuracy in the buffer packing routine.

* ``Test_N_VMBufUnpack``: Tests for accuracy in the buffer unpacking routine.
