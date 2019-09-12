#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_POWER_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_POWER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>

namespace stan {
namespace math {
/**
 * Returns the nth power of the specific matrix.
 *
 * @param M A square matrix.
 * @param n Exponent.
 * @return nth power of M. M^n = M * ... * M.
 * @throw std::invalid_argument if the exponent is negative or the matrix is not
 * square.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_power(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, const int n) {
  check_square("matrix_power", "M", M);
  if (n < 0)
    invalid_argument("matrix_power", "n", n, "is ", ", but must be >= 0!");
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result
      = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(M.rows(),
                                                                   M.rows());
  int nn = n;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MM = M;

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

}  // namespace math
}  // namespace stan
#endif
