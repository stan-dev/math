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
 * @tparam Ta scalar type matrix A
 * @tparam Tb scalar type matrix B
 * @tparam Cb Columns matrix B
 *
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return exponential of A multiplies B
 */
template <typename Ta, typename Tb, int Cb>
inline Eigen::Matrix<typename stan::return_type<Ta, Tb>::type, -1, Cb>
matrix_exp_multiply(const Eigen::Matrix<Ta, -1, -1>& A,
                    const Eigen::Matrix<Tb, -1, Cb>& B) {
  check_square("matrix_exp_multiply", "input matrix", A);
  if (A.size() == 0 && B.rows() == 0) {
    return Eigen::Matrix<typename stan::return_type_t<Ta, Tb>, -1, Cb>(
        0, B.cols());
  }

  check_multiplicable("matrix_exp_multiply", "A", A, "B", B);

  return multiply(matrix_exp(A), B);
}

}  // namespace math
}  // namespace stan

#endif
