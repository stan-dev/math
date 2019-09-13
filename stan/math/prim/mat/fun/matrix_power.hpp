#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_POWER_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_POWER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_finite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>

namespace stan {
namespace math {
/**
 * Returns the nth power of the specific matrix. M^n = M * M * ... * M.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows in matrix
 * @tparam C number of columns in matrix
 * @param[in] M a square matrix
 * @param[in] n exponent
 * @return nth power of M
 * @throw std::domain_error if the matrix contains NaNs or infinities.
 * @throw std::invalid_argument if the exponent is negative or the matrix is not
 * square.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> matrix_power(const Eigen::Matrix<T, R, C> &M,
                                           const int n) {
  check_square("matrix_power", "M", M);
  if (n < 0)
    invalid_argument("matrix_power", "n", n, "is ", ", but must be >= 0!");
  if (M.rows() == 0)
    invalid_argument("matrix_power", "M.rows()", M.rows(), "is ",
                     ", but must be > 0!");
  check_finite("matrix_power", "M", M);
  Eigen::Matrix<T, R, C> result
      = Eigen::Matrix<T, R, C>::Identity(M.rows(), M.cols());
  int nn = n;
  Eigen::Matrix<T, R, C> MM = M;
  while (nn > 0) {
    if (nn % 2 == 1) {
      result = result * MM;
      --nn;
    }
    MM = MM * MM;
    nn /= 2;
  }
  return result;
}

template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> operator^(const Eigen::Matrix<T, R, C> &M,
                                        const int n) {
  return matrix_power(M, n);
}

}  // namespace math
}  // namespace stan
#endif
