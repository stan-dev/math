// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2007-2011 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_CORE_H
#define EIGEN_CORE_H

// first thing Eigen does: stop the compiler from committing suicide
#include "Eigen/src/Core/util/DisableStupidWarnings.h"

#if defined(__CUDACC__) && !defined(EIGEN_NO_CUDA)
  #define EIGEN_CUDACC __CUDACC__
#endif

#if defined(__CUDA_ARCH__) && !defined(EIGEN_NO_CUDA)
  #define EIGEN_CUDA_ARCH __CUDA_ARCH__
#endif

#if defined(__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ >= 9)
#define EIGEN_CUDACC_VER  ((__CUDACC_VER_MAJOR__ * 10000) + (__CUDACC_VER_MINOR__ * 100))
#elif defined(__CUDACC_VER__)
#define EIGEN_CUDACC_VER __CUDACC_VER__
#else
#define EIGEN_CUDACC_VER 0
#endif

// Handle NVCC/CUDA/SYCL
#if defined(__CUDACC__) || defined(__SYCL_DEVICE_ONLY__)
  // Do not try asserts on CUDA and SYCL!
  #ifndef EIGEN_NO_DEBUG
  #define EIGEN_NO_DEBUG
  #endif

  #ifdef EIGEN_INTERNAL_DEBUGGING
  #undef EIGEN_INTERNAL_DEBUGGING
  #endif

  #ifdef EIGEN_EXCEPTIONS
  #undef EIGEN_EXCEPTIONS
  #endif

  // All functions callable from CUDA code must be qualified with __device__
  #ifdef __CUDACC__
    // Do not try to vectorize on CUDA and SYCL!
    #ifndef EIGEN_DONT_VECTORIZE
    #define EIGEN_DONT_VECTORIZE
    #endif

    #define EIGEN_DEVICE_FUNC __host__ __device__
    // We need cuda_runtime.h to ensure that that EIGEN_USING_STD_MATH macro
    // works properly on the device side
    #include <cuda_runtime.h>
  #else
    #define EIGEN_DEVICE_FUNC
  #endif

#else
  #define EIGEN_DEVICE_FUNC

#endif

// When compiling CUDA device code with NVCC, pull in math functions from the
// global namespace.  In host mode, and when device doee with clang, use the
// std versions.
#if defined(__CUDA_ARCH__) && defined(__NVCC__)
  #define EIGEN_USING_STD_MATH(FUNC) using ::FUNC;
#else
  #define EIGEN_USING_STD_MATH(FUNC) using std::FUNC;
#endif

#if (defined(_CPPUNWIND) || defined(__EXCEPTIONS)) && !defined(__CUDA_ARCH__) && !defined(EIGEN_EXCEPTIONS) && !defined(EIGEN_USE_SYCL)
  #define EIGEN_EXCEPTIONS
#endif

#ifdef EIGEN_EXCEPTIONS
  #include <new>
#endif

// then include this file where all our macros are defined. It's really important to do it first because
// it's where we do all the alignment settings (platform detection and honoring the user's will if he
// defined e.g. EIGEN_DONT_ALIGN) so it needs to be done before we do anything with vectorization.
#include "Eigen/src/Core/util/Macros.h"

// Disable the ipa-cp-clone optimization flag with MinGW 6.x or newer (enabled by default with -O3)
// See http://eigen.tuxfamily.org/bz/show_bug.cgi?id=556 for details.
#if EIGEN_COMP_MINGW && EIGEN_GNUC_AT_LEAST(4,6)
  #pragma GCC optimize ("-fno-ipa-cp-clone")
#endif

#include <complex>

// this include file manages BLAS and MKL related macros
// and inclusion of their respective header files
#include "Eigen/src/Core/util/MKL_support.h"

// if alignment is disabled, then disable vectorization. Note: EIGEN_MAX_ALIGN_BYTES is the proper check, it takes into
// account both the user's will (EIGEN_MAX_ALIGN_BYTES,EIGEN_DONT_ALIGN) and our own platform checks
#if EIGEN_MAX_ALIGN_BYTES==0
  #ifndef EIGEN_DONT_VECTORIZE
    #define EIGEN_DONT_VECTORIZE
  #endif
#endif

#if EIGEN_COMP_MSVC
  #include <malloc.h> // for _aligned_malloc -- need it regardless of whether vectorization is enabled
  #if (EIGEN_COMP_MSVC >= 1500) // 2008 or later
    // Remember that usage of defined() in a #define is undefined by the standard.
    // a user reported that in 64-bit mode, MSVC doesn't care to define _M_IX86_FP.
    #if (defined(_M_IX86_FP) && (_M_IX86_FP >= 2)) || EIGEN_ARCH_x86_64
      #define EIGEN_SSE2_ON_MSVC_2008_OR_LATER
    #endif
  #endif
#else
  // Remember that usage of defined() in a #define is undefined by the standard
  #if (defined __SSE2__) && ( (!EIGEN_COMP_GNUC) || EIGEN_COMP_ICC || EIGEN_GNUC_AT_LEAST(4,2) )
    #define EIGEN_SSE2_ON_NON_MSVC_BUT_NOT_OLD_GCC
  #endif
#endif

#ifndef EIGEN_DONT_VECTORIZE

  #if defined (EIGEN_SSE2_ON_NON_MSVC_BUT_NOT_OLD_GCC) || defined(EIGEN_SSE2_ON_MSVC_2008_OR_LATER)

    // Defines symbols for compile-time detection of which instructions are
    // used.
    // EIGEN_VECTORIZE_YY is defined if and only if the instruction set YY is used
    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_SSE
    #define EIGEN_VECTORIZE_SSE2

    // Detect sse3/ssse3/sse4:
    // gcc and icc defines __SSE3__, ...
    // there is no way to know about this on msvc. You can define EIGEN_VECTORIZE_SSE* if you
    // want to force the use of those instructions with msvc.
    #ifdef __SSE3__
      #define EIGEN_VECTORIZE_SSE3
    #endif
    #ifdef __SSSE3__
      #define EIGEN_VECTORIZE_SSSE3
    #endif
    #ifdef __SSE4_1__
      #define EIGEN_VECTORIZE_SSE4_1
    #endif
    #ifdef __SSE4_2__
      #define EIGEN_VECTORIZE_SSE4_2
    #endif
    #ifdef __AVX__
      #define EIGEN_VECTORIZE_AVX
      #define EIGEN_VECTORIZE_SSE3
      #define EIGEN_VECTORIZE_SSSE3
      #define EIGEN_VECTORIZE_SSE4_1
      #define EIGEN_VECTORIZE_SSE4_2
    #endif
    #ifdef __AVX2__
      #define EIGEN_VECTORIZE_AVX2
    #endif
    #ifdef __FMA__
      #define EIGEN_VECTORIZE_FMA
    #endif
    #if defined(__AVX512F__) && defined(EIGEN_ENABLE_AVX512)
      #define EIGEN_VECTORIZE_AVX512
      #define EIGEN_VECTORIZE_AVX2
      #define EIGEN_VECTORIZE_AVX
      #define EIGEN_VECTORIZE_FMA
      #ifdef __AVX512DQ__
        #define EIGEN_VECTORIZE_AVX512DQ
      #endif
      #ifdef __AVX512ER__
        #define EIGEN_VECTORIZE_AVX512ER
      #endif
    #endif

    // include files

    // This extern "C" works around a MINGW-w64 compilation issue
    // https://sourceforge.net/tracker/index.php?func=detail&aid=3018394&group_id=202880&atid=983354
    // In essence, intrin.h is included by windows.h and also declares intrinsics (just as emmintrin.h etc. below do).
    // However, intrin.h uses an extern "C" declaration, and g++ thus complains of duplicate declarations
    // with conflicting linkage.  The linkage for intrinsics doesn't matter, but at that stage the compiler doesn't know;
    // so, to avoid compile errors when windows.h is included after Eigen/Core, ensure intrinsics are extern "C" here too.
    // notice that since these are C headers, the extern "C" is theoretically needed anyways.
    extern "C" {
      // In theory we should only include immintrin.h and not the other *mmintrin.h header files directly.
      // Doing so triggers some issues with ICC. However old gcc versions seems to not have this file, thus:
      #if EIGEN_COMP_ICC >= 1110
        #include <immintrin.h>
      #else
        #include <mmintrin.h>
        #include <emmintrin.h>
        #include <xmmintrin.h>
        #ifdef  EIGEN_VECTORIZE_SSE3
        #include <pmmintrin.h>
        #endif
        #ifdef EIGEN_VECTORIZE_SSSE3
        #include <tmmintrin.h>
        #endif
        #ifdef EIGEN_VECTORIZE_SSE4_1
        #include <smmintrin.h>
        #endif
        #ifdef EIGEN_VECTORIZE_SSE4_2
        #include <nmmintrin.h>
        #endif
        #if defined(EIGEN_VECTORIZE_AVX) || defined(EIGEN_VECTORIZE_AVX512)
        #include <immintrin.h>
        #endif
      #endif
    } // end extern "C"
  #elif defined __VSX__
    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_VSX
    #include <altivec.h>
    // We need to #undef all these ugly tokens defined in <altivec.h>
    // => use __vector instead of vector
    #undef bool
    #undef vector
    #undef pixel
  #elif defined __ALTIVEC__
    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_ALTIVEC
    #include <altivec.h>
    // We need to #undef all these ugly tokens defined in <altivec.h>
    // => use __vector instead of vector
    #undef bool
    #undef vector
    #undef pixel
  #elif (defined  __ARM_NEON) || (defined __ARM_NEON__)
    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_NEON
    #include <arm_neon.h>
  #elif (defined __s390x__ && defined __VEC__)
    #define EIGEN_VECTORIZE
    #define EIGEN_VECTORIZE_ZVECTOR
    #include <vecintrin.h>
  #endif
#endif

#if defined(__F16C__) && !defined(EIGEN_COMP_CLANG)
  // We can use the optimized fp16 to float and float to fp16 conversion routines
  #define EIGEN_HAS_FP16_C
#endif

#if defined __CUDACC__
  #define EIGEN_VECTORIZE_CUDA
  #include <vector_types.h>
  #if EIGEN_CUDACC_VER >= 70500
    #define EIGEN_HAS_CUDA_FP16
  #endif
#endif

#if defined EIGEN_HAS_CUDA_FP16
  #include <host_defines.h>
  #include <cuda_fp16.h>
#endif

#if (defined _OPENMP) && (!defined EIGEN_DONT_PARALLELIZE)
  #define EIGEN_HAS_OPENMP
#endif

#ifdef EIGEN_HAS_OPENMP
#include <omp.h>
#endif

// MSVC for windows mobile does not have the errno.h file
#if !(EIGEN_COMP_MSVC && EIGEN_OS_WINCE) && !EIGEN_COMP_ARM
#define EIGEN_HAS_ERRNO
#endif

#ifdef EIGEN_HAS_ERRNO
#include <cerrno>
#endif
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <functional>
#include <sstream>
#ifndef EIGEN_NO_IO
  #include <iosfwd>
#endif
#include <cstring>
#include <string>
#include <limits>
#include <climits> // for CHAR_BIT
// for min/max:
#include <algorithm>

// for std::is_nothrow_move_assignable
#ifdef EIGEN_INCLUDE_TYPE_TRAITS
#include <type_traits>
#endif

// for outputting debug info
#ifdef EIGEN_DEBUG_ASSIGN
#include <iostream>
#endif

// required for __cpuid, needs to be included after cmath
#if EIGEN_COMP_MSVC && EIGEN_ARCH_i386_OR_x86_64 && !EIGEN_OS_WINCE
  #include <intrin.h>
#endif

/** \brief Namespace containing all symbols from the %Eigen library. */
namespace Eigen {

inline static const char *SimdInstructionSetsInUse(void) {
#if defined(EIGEN_VECTORIZE_AVX512)
  return "AVX512, FMA, AVX2, AVX, SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2";
#elif defined(EIGEN_VECTORIZE_AVX)
  return "AVX SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2";
#elif defined(EIGEN_VECTORIZE_SSE4_2)
  return "SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2";
#elif defined(EIGEN_VECTORIZE_SSE4_1)
  return "SSE, SSE2, SSE3, SSSE3, SSE4.1";
#elif defined(EIGEN_VECTORIZE_SSSE3)
  return "SSE, SSE2, SSE3, SSSE3";
#elif defined(EIGEN_VECTORIZE_SSE3)
  return "SSE, SSE2, SSE3";
#elif defined(EIGEN_VECTORIZE_SSE2)
  return "SSE, SSE2";
#elif defined(EIGEN_VECTORIZE_ALTIVEC)
  return "AltiVec";
#elif defined(EIGEN_VECTORIZE_VSX)
  return "VSX";
#elif defined(EIGEN_VECTORIZE_NEON)
  return "ARM NEON";
#elif defined(EIGEN_VECTORIZE_ZVECTOR)
  return "S390X ZVECTOR";
#else
  return "None";
#endif
}

} // end namespace Eigen

#if defined EIGEN2_SUPPORT_STAGE40_FULL_EIGEN3_STRICTNESS || defined EIGEN2_SUPPORT_STAGE30_FULL_EIGEN3_API || defined EIGEN2_SUPPORT_STAGE20_RESOLVE_API_CONFLICTS || defined EIGEN2_SUPPORT_STAGE10_FULL_EIGEN2_API || defined EIGEN2_SUPPORT
// This will generate an error message:
#error Eigen2-support is only available up to version 3.2. Please go to "http://eigen.tuxfamily.org/index.php?title=Eigen2" for further information
#endif

namespace Eigen {

// we use size_t frequently and we'll never remember to prepend it with std:: everytime just to
// ensure QNX/QCC support
using std::size_t;
// gcc 4.6.0 wants std:: for ptrdiff_t
using std::ptrdiff_t;

}

/** \defgroup Core_Module Core module
  * This is the main module of Eigen providing dense matrix and vector support
  * (both fixed and dynamic size) with all the features corresponding to a BLAS library
  * and much more...
  *
  * \code
  * #include <Eigen/Core>
  * \endcode
  */

#include "Eigen/src/Core/util/Constants.h"
#include "Eigen/src/Core/util/Meta.h"
#include "stan/math/prim/plugins/ForwardDeclarations.h"
#include "Eigen/src/Core/util/StaticAssert.h"
#include "Eigen/src/Core/util/XprHelper.h"
#include "Eigen/src/Core/util/Memory.h"

#include "Eigen/src/Core/NumTraits.h"
#include "Eigen/src/Core/MathFunctions.h"
#include "Eigen/src/Core/GenericPacketMath.h"
#include "Eigen/src/Core/MathFunctionsImpl.h"
#include "Eigen/src/Core/arch/Default/ConjHelper.h"

#if defined EIGEN_VECTORIZE_AVX512
  #include "Eigen/src/Core/arch/SSE/PacketMath.h"
  #include "Eigen/src/Core/arch/SSE/MathFunctions.h"
  #include "Eigen/src/Core/arch/AVX/PacketMath.h"
  #include "Eigen/src/Core/arch/AVX/MathFunctions.h"
  #include "Eigen/src/Core/arch/AVX512/PacketMath.h"
  #include "Eigen/src/Core/arch/AVX512/MathFunctions.h"
#elif defined EIGEN_VECTORIZE_AVX
  // Use AVX for floats and doubles, SSE for integers
  #include "Eigen/src/Core/arch/SSE/PacketMath.h"
  #include "Eigen/src/Core/arch/SSE/Complex.h"
  #include "Eigen/src/Core/arch/SSE/MathFunctions.h"
  #include "Eigen/src/Core/arch/AVX/PacketMath.h"
  #include "Eigen/src/Core/arch/AVX/MathFunctions.h"
  #include "Eigen/src/Core/arch/AVX/Complex.h"
  #include "Eigen/src/Core/arch/AVX/TypeCasting.h"
  #include "Eigen/src/Core/arch/SSE/TypeCasting.h"
#elif defined EIGEN_VECTORIZE_SSE
  #include "Eigen/src/Core/arch/SSE/PacketMath.h"
  #include "Eigen/src/Core/arch/SSE/MathFunctions.h"
  #include "Eigen/src/Core/arch/SSE/Complex.h"
  #include "Eigen/src/Core/arch/SSE/TypeCasting.h"
#elif defined(EIGEN_VECTORIZE_ALTIVEC) || defined(EIGEN_VECTORIZE_VSX)
  #include "Eigen/src/Core/arch/AltiVec/PacketMath.h"
  #include "Eigen/src/Core/arch/AltiVec/MathFunctions.h"
  #include "Eigen/src/Core/arch/AltiVec/Complex.h"
#elif defined EIGEN_VECTORIZE_NEON
  #include "Eigen/src/Core/arch/NEON/PacketMath.h"
  #include "Eigen/src/Core/arch/NEON/MathFunctions.h"
  #include "Eigen/src/Core/arch/NEON/Complex.h"
#elif defined EIGEN_VECTORIZE_ZVECTOR
  #include "Eigen/src/Core/arch/ZVector/PacketMath.h"
  #include "Eigen/src/Core/arch/ZVector/MathFunctions.h"
  #include "Eigen/src/Core/arch/ZVector/Complex.h"
#endif

// Half float support
#include "Eigen/src/Core/arch/CUDA/Half.h"
#include "Eigen/src/Core/arch/CUDA/PacketMathHalf.h"
#include "Eigen/src/Core/arch/CUDA/TypeCasting.h"

#if defined EIGEN_VECTORIZE_CUDA
  #include "Eigen/src/Core/arch/CUDA/PacketMath.h"
  #include "Eigen/src/Core/arch/CUDA/MathFunctions.h"
#endif

#include "Eigen/src/Core/arch/Default/Settings.h"

#include "Eigen/src/Core/functors/TernaryFunctors.h"
#include "Eigen/src/Core/functors/BinaryFunctors.h"
#include "Eigen/src/Core/functors/UnaryFunctors.h"
#include "Eigen/src/Core/functors/NullaryFunctors.h"
#include "Eigen/src/Core/functors/StlFunctors.h"
#include "Eigen/src/Core/functors/AssignmentFunctors.h"

// Specialized functors to enable the processing of complex numbers
// on CUDA devices
#include "Eigen/src/Core/arch/CUDA/Complex.h"

#include "Eigen/src/Core/IO.h"
#include "Eigen/src/Core/DenseCoeffsBase.h"
#include "Eigen/src/Core/DenseBase.h"
#include "Eigen/src/Core/MatrixBase.h"
#include "Eigen/src/Core/EigenBase.h"

#include "Eigen/src/Core/Product.h"
#include "stan/math/prim/plugins/CoreEvaluators.h"
#include "Eigen/src/Core/AssignEvaluator.h"

#ifndef EIGEN_PARSED_BY_DOXYGEN // work around Doxygen bug triggered by Assign.h r814874
                                // at least confirmed with Doxygen 1.5.5 and 1.5.6
  #include "Eigen/src/Core/Assign.h"
#endif

#include "Eigen/src/Core/ArrayBase.h"
#include "Eigen/src/Core/util/BlasUtil.h"
#include "Eigen/src/Core/DenseStorage.h"
#include "Eigen/src/Core/NestByValue.h"

// #include "Eigen/src/Core/ForceAlignedAccess.h"

#include "Eigen/src/Core/ReturnByValue.h"
#include "Eigen/src/Core/NoAlias.h"
#include "Eigen/src/Core/PlainObjectBase.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/Array.h"
#include "Eigen/src/Core/CwiseTernaryOp.h"
#include "Eigen/src/Core/CwiseBinaryOp.h"
#include "Eigen/src/Core/CwiseUnaryOp.h"
#include "Eigen/src/Core/CwiseNullaryOp.h"
#include "stan/math/prim/plugins/CwiseUnaryView.h"
#include "Eigen/src/Core/SelfCwiseBinaryOp.h"
#include "Eigen/src/Core/Dot.h"
#include "Eigen/src/Core/StableNorm.h"
#include "Eigen/src/Core/Stride.h"
#include "Eigen/src/Core/MapBase.h"
#include "Eigen/src/Core/Map.h"
#include "Eigen/src/Core/Ref.h"
#include "Eigen/src/Core/Block.h"
#include "Eigen/src/Core/VectorBlock.h"
#include "Eigen/src/Core/Transpose.h"
#include "Eigen/src/Core/DiagonalMatrix.h"
#include "Eigen/src/Core/Diagonal.h"
#include "Eigen/src/Core/DiagonalProduct.h"
#include "Eigen/src/Core/Redux.h"
#include "Eigen/src/Core/Visitor.h"
#include "Eigen/src/Core/Fuzzy.h"
#include "Eigen/src/Core/Swap.h"
#include "Eigen/src/Core/CommaInitializer.h"
#include "Eigen/src/Core/GeneralProduct.h"
#include "Eigen/src/Core/Solve.h"
#include "Eigen/src/Core/Inverse.h"
#include "Eigen/src/Core/SolverBase.h"
#include "Eigen/src/Core/PermutationMatrix.h"
#include "Eigen/src/Core/Transpositions.h"
#include "Eigen/src/Core/TriangularMatrix.h"
#include "Eigen/src/Core/SelfAdjointView.h"
#include "Eigen/src/Core/products/GeneralBlockPanelKernel.h"
#include "Eigen/src/Core/products/Parallelizer.h"
#include "Eigen/src/Core/ProductEvaluators.h"
#include "Eigen/src/Core/products/GeneralMatrixVector.h"
#include "Eigen/src/Core/products/GeneralMatrixMatrix.h"
#include "Eigen/src/Core/SolveTriangular.h"
#include "Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h"
#include "Eigen/src/Core/products/SelfadjointMatrixVector.h"
#include "Eigen/src/Core/products/SelfadjointMatrixMatrix.h"
#include "Eigen/src/Core/products/SelfadjointProduct.h"
#include "Eigen/src/Core/products/SelfadjointRank2Update.h"
#include "Eigen/src/Core/products/TriangularMatrixVector.h"
#include "Eigen/src/Core/products/TriangularMatrixMatrix.h"
#include "Eigen/src/Core/products/TriangularSolverMatrix.h"
#include "Eigen/src/Core/products/TriangularSolverVector.h"
#include "Eigen/src/Core/BandMatrix.h"
#include "Eigen/src/Core/CoreIterators.h"
#include "Eigen/src/Core/ConditionEstimator.h"

#include "Eigen/src/Core/BooleanRedux.h"
#include "Eigen/src/Core/Select.h"
#include "Eigen/src/Core/VectorwiseOp.h"
#include "Eigen/src/Core/Random.h"
#include "Eigen/src/Core/Replicate.h"
#include "Eigen/src/Core/Reverse.h"
#include "Eigen/src/Core/ArrayWrapper.h"

#ifdef EIGEN_USE_BLAS
#include "Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h"
#include "Eigen/src/Core/products/GeneralMatrixVector_BLAS.h"
#include "Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h"
#include "Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h"
#include "Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h"
#include "Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h"
#include "Eigen/src/Core/products/TriangularMatrixVector_BLAS.h"
#include "Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h"
#endif // EIGEN_USE_BLAS

#ifdef EIGEN_USE_MKL_VML
#include "Eigen/src/Core/Assign_MKL.h"
#endif

#include "Eigen/src/Core/GlobalFunctions.h"

#include "Eigen/src/Core/util/ReenableStupidWarnings.h"

#endif // EIGEN_CORE_H
