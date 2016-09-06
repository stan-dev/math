#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_IS_SYMMETRIC_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_IS_SYMMETRIC_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/constraint_tolerance.hpp>
#include <stan/math/prim/mat/meta/index_type.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>


namespace stan {
    namespace math {
        
        /**
         * Return a boolean for whether a matrix is symmetric
         * or not.
         * @param A A matrix
         * @return Boolean: True is the matrix is symmetric, false
         * otherwise.
         */
        
        template <typename T>
        inline
        bool
        is_symmetric(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A) {
            
            typedef typename index_type<Matrix<T, Dynamic, Dynamic> >::type
            size_type;
            
            size_type k = A.rows();
            if (k == 1) return true;
            for (size_type m = 0; m < k; ++m) {
                for (size_type n = m + 1; n < k; ++n) {
                    if (!(fabs(value_of(A(m, n)) - value_of(A(n, m)))
                          <= CONSTRAINT_TOLERANCE)) return false;
                }
            }
            return true;
        }
        
    } // end namespace math
} // namespace stan

#endif