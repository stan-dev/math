#ifndef STAN_MATH_PRIM_MAT_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_PRIM_MAT_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/rev/mat/fun/matrix_exp_action_handler.hpp>

namespace stan {
namespace math {

/**
 * Return product of exp(At) and B, where A is a NxN double matrix,
 * B is a NxCb double matrix, and t is a double
 * @tparam N Rows and cols matrix A, also rows of matrix B
 * @tparam Cb Columns matrix B
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplies B
 */
template <int N, int Cb>
inline Eigen::Matrix<double, N, Cb> scale_matrix_exp_multiply(
    const double& t, const Eigen::Matrix<double, N, N>& A,
    const Eigen::Matrix<double, N, Cb>& B) {
  Eigen::Matrix<double, N, Cb> expAB;
  stan::math::matrix_exp_action_handler handle;
  expAB = handle.action(A, B, t);
  return expAB;
}

}  // namespace math
}  // namespace stan
#endif
