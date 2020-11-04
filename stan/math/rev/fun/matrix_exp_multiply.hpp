#ifndef STAN_MATH_REV_FUN_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/matrix_exp.hpp>

namespace stan {
namespace math {

/**
 * Wrapper of matrix_exp_action function for a more literal name
 *
 * @tparam Ta type of the matrix A
 * @tparam Tb type of the matrix B
 *
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return exponential of A multiplies B
 */
template <typename Ta, typename Tb, require_all_eigen_t<Ta, Tb>* = nullptr,
          require_any_st_autodiff<Ta, Tb>* = nullptr>
inline Eigen::Matrix<return_type_t<Ta, Tb>, -1, Tb::ColsAtCompileTime>
matrix_exp_multiply(const Ta& A, const Tb& B) {
  check_square("matrix_exp_multiply", "input matrix", A);
  check_multiplicable("matrix_exp_multiply", "A", A, "B", B);
  if (A.size() == 0) {
    return {0, B.cols()};
  }

  return multiply(matrix_exp(A), B);
}

}  // namespace math
}  // namespace stan

#endif
