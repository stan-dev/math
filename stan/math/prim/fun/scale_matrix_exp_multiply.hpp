#ifndef STAN_MATH_PRIM_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/matrix_exp.hpp>
#include <stan/math/prim/fun/matrix_exp_action_handler.hpp>
#include <stan/math/prim/fun/multiply.hpp>

namespace stan {
namespace math {

/**
 * Return product of exp(At) and B, where A is a NxN double matrix,
 * B is a NxCb double matrix, and t is a double.
 *
 * Specialized for double values for efficiency.
 *
 * @tparam EigMat1 type of the first matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_vt<std::is_arithmetic, EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<double, Eigen::Dynamic, EigMat2::ColsAtCompileTime>
scale_matrix_exp_multiply(const double& t, const EigMat1& A, const EigMat2& B) {
  check_square("scale_matrix_exp_multiply", "input matrix", A);
  check_multiplicable("scale_matrix_exp_multiply", "A", A, "B", B);
  if (A.size() == 0) {
    return {0, B.cols()};
  }

  return matrix_exp_action_handler().action(A, B, t);
}

/**
 * Return product of exp(At) and B, where A is a NxN matrix,
 * B is a NxCb matrix and t is a scalar.
 *
 * Generic implementation when arguments are not all double.
 *
 * @tparam Tt type of \c t
 * @tparam EigMat1 type of the first matrix
 * @tparam EigMat2 type of the second matrix
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplied by B
 */
template <typename Tt, typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_any_autodiff_t<Tt, value_type_t<EigMat1>,
                                 value_type_t<EigMat2>>* = nullptr>
inline Eigen::Matrix<return_type_t<Tt, EigMat1, EigMat2>, Eigen::Dynamic,
                     EigMat2::ColsAtCompileTime>
scale_matrix_exp_multiply(const Tt& t, const EigMat1& A, const EigMat2& B) {
  check_square("scale_matrix_exp_multiply", "input matrix", A);
  check_multiplicable("scale_matrix_exp_multiply", "A", A, "B", B);
  if (A.size() == 0) {
    return {0, B.cols()};
  }

  return multiply(matrix_exp(multiply(A, t)), B);
}

}  // namespace math
}  // namespace stan
#endif
