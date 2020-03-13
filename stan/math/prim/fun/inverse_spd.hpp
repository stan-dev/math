#ifndef STAN_MATH_PRIM_FUN_INVERSE_SPD_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_SPD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified symmetric, pos/neg-definite matrix.
 *
 * @tparam T type of elements in the matrix
 * @param m specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not symmetric.
 * @throw std::domain_error if the matrix is not positive definite.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inverse_spd(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  using Eigen::Dynamic;
  using Eigen::LDLT;
  using Eigen::Matrix;
  check_symmetric("inverse_spd", "m", m);
  if (m.size() == 0) {
    return {};
  }

  Matrix<T, Dynamic, Dynamic> mmt = T(0.5) * (m + m.transpose());
  LDLT<Matrix<T, Dynamic, Dynamic> > ldlt(mmt);
  if (ldlt.info() != Eigen::Success) {
    throw_domain_error("invese_spd", "LDLT factor failed", "", "");
  }
  if (!ldlt.isPositive()) {
    throw_domain_error("invese_spd", "matrix not positive definite", "", "");
  }
  Matrix<T, Dynamic, 1> diag_ldlt = ldlt.vectorD();
  for (int i = 0; i < diag_ldlt.size(); ++i) {
    if (diag_ldlt(i) <= 0) {
      throw_domain_error("invese_spd", "matrix not positive definite", "", "");
    }
  }
  return ldlt.solve(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(
      m.rows(), m.cols()));
}

}  // namespace math
}  // namespace stan

#endif
