#ifndef STAN_MATH_PRIM_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/matrix_exp_action_handler.hpp>

namespace stan {
namespace math {

/**
 * Return product of exp(At) and B, where A is a NxN double matrix,
 * B is a NxCb double matrix, and t is a double.
 *
 * Specialized for double values for efficiency.
 *
 * @tparam Cb number of columns in matrix B, can be Eigen::Dynamic
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplies B
 */
template <int Cb>
inline Eigen::Matrix<double, -1, Cb> scale_matrix_exp_multiply(
    const double& t, const Eigen::MatrixXd& A,
    const Eigen::Matrix<double, -1, Cb>& B) {
  check_square("scale_matrix_exp_multiply", "input matrix", A);
  if (A.size() == 0 && B.rows() == 0) {
    return Eigen::Matrix<double, -1, Cb>(0, B.cols());
  }

  check_multiplicable("scale_matrix_exp_multiply", "A", A, "B", B);

  return matrix_exp_action_handler().action(A, B, t);
}

/**
 * Return product of exp(At) and B, where A is a NxN matrix,
 * B is a NxCb matrix and t is a scalar.
 *
 * Generic implementation when arguments are not double.
 *
 * @tparam Ta scalar type matrix A
 * @tparam Tb scalar type matrix B
 * @tparam Cb number of columns in matrix B, can be Eigen::Dynamic
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplies B
 */
template <typename Tt, typename Ta, typename Tb, int Cb>
inline Eigen::Matrix<return_type_t<Tt, Ta, Tb>, -1, Cb>
scale_matrix_exp_multiply(const Tt& t, const Eigen::Matrix<Ta, -1, -1>& A,
                          const Eigen::Matrix<Tb, -1, Cb>& B) {
  check_square("scale_matrix_exp_multiply", "input matrix", A);
  if (A.size() == 0 && B.rows() == 0) {
    return Eigen::Matrix<return_type_t<Tt, Ta, Tb>, -1, Cb>(0, B.cols());
  }

  check_multiplicable("scale_matrix_exp_multiply", "A", A, "B", B);

  return multiply(matrix_exp(multiply(A, t)), B);
}

}  // namespace math
}  // namespace stan
#endif
