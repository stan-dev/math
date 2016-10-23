#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_PADE_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_PADE_HPP

#include <stan/math/prim/mat/fun/MatrixExponential.h>

namespace stan {
  namespace math {

        /**
         * Computes the matrix exponential, using a Pade
         * approximation, coupled with scaling and
         * squaring.
         *
         * @tparam MatrixType scalar type of the elements
         * in the input matrix.
         * @param[in] arg matrix to exponentiate.
         * @param[out] Matrix exponential of input matrix.
         */
        template <typename MatrixType>
        MatrixType
        matrix_exp_pade(const MatrixType& arg) {
          MatrixType result;

        #if LDBL_MANT_DIG > 112  // rarely happens
          typedef typename Eigen::traits<MatrixType>::Scalar Scalar;
          typedef typename Eigen::NumTraits<Scalar>::Real RealScalar;
          typedef typename std::complex<RealScalar> ComplexScalar;
          if (sizeof(RealScalar) > 14) {
            result = arg.matrixFunction(
              internal::stem_function_exp<ComplexScalar>);
            return;
            }
        #endif
          MatrixType U, V;
          int squarings;
          Eigen:: matrix_exp_computeUV<MatrixType>::run(
            arg, U, V, squarings, arg(0, 0));  // Pade approximant is
                                               // (U+V) / (-U+V)
          MatrixType numer = U + V;
          MatrixType denom = -U + V;
          result = denom.partialPivLu().solve(numer);
          for (int i = 0; i < squarings; i++)
            result *= result;  // undo scaling by repeated squaring

          return result;
        }

    }
}

#endif
