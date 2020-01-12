#ifndef STAN_MATH_PRIM_FUN_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/matrix_exp_action_handler.hpp>

namespace stan {
namespace math {

/**
 * Return product of exp(A) and B, where A is a NxN double matrix,
 * B is a NxCb double matrix, and t is a double
 *
 * @tparam Cb number of columns in matrix B, can be Eigen::Dynamic
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return exponential of A multiplies B
 */
template <int Cb>
inline Eigen::Matrix<double, -1, Cb> matrix_exp_multiply(
    const Eigen::MatrixXd& A, const Eigen::Matrix<double, -1, Cb>& B) {
  check_square("matrix_exp_multiply", "input matrix", A);
  if (A.size() == 0 && B.rows() == 0) {
    return Eigen::Matrix<double, -1, Cb>(0, B.cols());
  }

  check_multiplicable("matrix_exp_multiply", "A", A, "B", B);

  return matrix_exp_action_handler().action(A, B);
}

}  // namespace math
}  // namespace stan

#endif
