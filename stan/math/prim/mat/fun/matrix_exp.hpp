#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <unsupported/Eigen/MatrixFunctions>

namespace stan {
    namespace math {
        
        /**
         * Return the matrix exponential.
         * @param A A matrix
         * @return Matrix exponential. 
         */
        
        template <typename T>
        inline
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
        matrix_exp(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A) {
            return A.exp();
        }
        
    }
}
#endif