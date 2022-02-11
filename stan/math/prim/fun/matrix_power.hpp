#ifndef STAN_MATH_PRIM_FUN_MATRIX_POWER_HPP
#define STAN_MATH_PRIM_FUN_MATRIX_POWER_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the nth power of the specific matrix. M^n = M * M * ... * M.
 *
 * @tparam T type of the matrix
 *
 * @param[in] M a square matrix
 * @param[in] n exponent
 * @return nth power of M
 * @throw std::domain_error if the matrix contains NaNs or infinities.
 * @throw std::invalid_argument if the exponent is negative or the matrix is not
 * square.
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
matrix_power(const EigMat& M, const int n) {
  using T = value_type_t<EigMat>;
  constexpr int R = EigMat::RowsAtCompileTime;
  constexpr int C = EigMat::ColsAtCompileTime;

  check_square("matrix_power", "M", M);
  check_nonnegative("matrix_power", "n", n);
  Eigen::Matrix<T, R, C> MM = M;
  check_finite("matrix_power", "M", MM);
  if (n == 0)
    return Eigen::Matrix<T, R, C>::Identity(M.rows(), M.cols());
  Eigen::Matrix<T, R, C> result = MM;
  for (int nn = n - 1; nn > 0; nn /= 2) {
    if (nn % 2 == 1) {
      result = result * MM;
      --nn;
    }
    MM = MM * MM;
  }
  return result;
}

template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
operator^(const EigMat& M, const int n) {
  return matrix_power(M, n);
}

}  // namespace math
}  // namespace stan

#endif
