#ifndef STAN_MATH_PRIM_FUN_INVERSE_SPD_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_SPD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified symmetric, pos/neg-definite matrix.
 *
 * @tparam EigMat type of elements in the matrix
 * @param m specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not symmetric.
 * @throw std::domain_error if the matrix is not positive definite.
 */
template <typename EigMat>
inline Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>
inverse_spd(const EigMat& m) {
  using Eigen::Dynamic;
  using Eigen::LDLT;
  using Eigen::Matrix;
  using Scalar = value_type_t<EigMat>;
  if (m.size() == 0) {
    return {};
  }
  const Eigen::Ref<const plain_type_t<EigMat>>& m_ref = m;
  check_symmetric("inverse_spd", "m", m_ref);
  plain_type_t<EigMat> mmt = 0.5 * (m_ref + m_ref.transpose());
  LDLT<plain_type_t<EigMat>> ldlt(mmt);
  if (ldlt.info() != Eigen::Success) {
    throw_domain_error("invese_spd", "LDLT factor failed", "", "");
  }
  if (!ldlt.isPositive()) {
    throw_domain_error("invese_spd", "matrix not positive definite", "", "");
  }
  Matrix<Scalar, Dynamic, 1> diag_ldlt = ldlt.vectorD();
  check_positive("inverse_spd", "matrix not positive definite", diag_ldlt);

  return ldlt.solve(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(
          m.rows(), m.cols()));
}

}  // namespace math
}  // namespace stan

#endif
