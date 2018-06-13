#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/prim/mat/fun/matrix_exp_action_handler.hpp>

namespace stan {
namespace math {

/**
 * Return product of exp(A) and B, where A is a NxN double matrix,
 * B is a NxCb double matrix, and t is a double
 * @tparam Cb Columns matrix B
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return exponential of A multiplies B
 */
template <int Cb>
inline Eigen::Matrix<double, -1, Cb> matrix_exp_multiply(
    const Eigen::MatrixXd& A,
    const Eigen::Matrix<double, -1, Cb>& B) {
  return matrix_exp_action_handler().action(A, B);
}

}  // namespace math
}  // namespace stan
#endif
