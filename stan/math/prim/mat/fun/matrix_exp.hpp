#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP

#include <stan/math/prim/mat/fun/matrix_exp_pade.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_2x2.hpp>

namespace stan {
  namespace math {

    /**
     * Return the matrix exponential of the input
     * matrix.
     *
     * @tparam T type of scalar of the elements of
     * input matrix.
     * @param[in] A Matrix to exponentiate.
     * @return Matrix exponential.
     */
    template <typename T>
    inline
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    matrix_exp(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A) {
      check_nonzero_size("matrix_exp", "input matrix", A);
      check_square("matrix_exp", "input matrix", A);

      return (A.cols() == 2
        && square(value_of(A(0, 0)) - value_of(A(1, 1)))
        + 4 * value_of(A(0, 1)) * value_of(A(1, 0)) > 0)
        ? matrix_exp_2x2(A) : matrix_exp_pade(A);
    }
  }
}
#endif
