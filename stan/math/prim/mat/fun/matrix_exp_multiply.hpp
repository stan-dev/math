#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/rev/mat/fun/matrix_exp_action_handler.hpp>

namespace stan {
namespace math {

/**
 * Return product of exp(A) and B, where A is a NxN double matrix,
 * B is a NxCb double matrix, and t is a double
 * @tparam N Rows and cols matrix A, also rows of matrix B
 * @tparam Cb Columns matrix B
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return exponential of A multiplies B
 */
  template <int N, int Cb>
  inline Eigen::Matrix<double, N, Cb>
  matrix_exp_multiply(const Eigen::Matrix<double, N, N>& A,
                      const Eigen::Matrix<double, N, Cb>& B) {
    Eigen::Matrix<double, N, Cb> expAB;
    stan::math::matrix_exp_action_handler handle;
    expAB = handle.action(A, B);
    return expAB;
  }

}  // namespace math
}  // namespace stan
#endif
