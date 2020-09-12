#ifndef STAN_MATH_PRIM_FUN_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
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
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_all_st_same<double, EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<double, Eigen::Dynamic, EigMat2::ColsAtCompileTime>
matrix_exp_multiply(const EigMat1& A, const EigMat2& B) {
  check_square("matrix_exp_multiply", "input matrix", A);
  check_multiplicable("matrix_exp_multiply", "A", A, "B", B);
  if (A.size() == 0) {
    return {0, B.cols()};
  }

  return matrix_exp_action_handler().action(A, B);
}

}  // namespace math
}  // namespace stan

#endif
