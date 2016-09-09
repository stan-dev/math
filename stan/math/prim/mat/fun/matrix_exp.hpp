#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP

#include <unsupported/Eigen/MatrixFunctions>
#include <stan/math/prim/mat/fun/MatrixExponential.h>
#include <stan/math/prim/mat/fun/matrix_exp_2x2.hpp>

namespace stan {
    namespace math {

        /**
         * Return the matrix exponential of the input
         * matrix. 
         *
         * @param[in] A A matrix
         * @return Matrix exponential. 
         */

        template <typename T>
        inline
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
        matrix_exp(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A) {
            check_nonzero_size("matrix_exp", "input matrix", A);
            check_square("matrix_exp", "input matrix", A);

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B;

            if (A.cols() == 2) B = matrix_exp_compute_2x2(A);
            else
               B = matrix_exp_compute(A);

            return B;
       }

    }
}
#endif
