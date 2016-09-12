#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_2X2_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_2X2_HPP

#include <unsupported/Eigen/MatrixFunctions>
#include <stan/math/prim/scal/err/check_size_match.hpp>

namespace stan {
    namespace math {

        /**
         * Return the matrix exponential of a 2x2 matrix. Reference for 
         * algorithm: http://mathworld.wolfram.com/MatrixExponential.html
         *
         * @param[in] A A 2x2 matrix
         * @return Matrix exponential of A.
         */
        template <typename T>
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
        matrix_exp_compute_2x2(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A) {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B;
            T a = A(0, 0), b = A(0, 1), c = A(1, 0), d = A(1, 1), delta;
            delta = sqrt((a - d) * (a - d) + 4 * b * c);

			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> B(2, 2);
			T half_delta = 0.5 * delta;
			T cosh_half_delta = cosh(half_delta);
			T sinh_half_delta = sinh(half_delta);
			T exp_half_a_plus_d = exp(0.5 * (a + d));
			
			B(0, 0) = exp_half_a_plus_d * (delta * cosh_half_delta + 
					  (a - d) * sinh_half_delta;
			B(0, 1) = 2 * b * exp_half_a_plus_d * sinh_half_delta;
			B(1, 0) = 2 * c * exp_half_a_plus_d * sinh_half_delta;
			B(1, 1) = exp_half_a_plus_d;

            return 1/delta * B;
        }

    }
}
#endif
