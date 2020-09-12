#ifndef STAN_MATH_FWD_FUN_INVERSE_HPP
#define STAN_MATH_FWD_FUN_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/inverse.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>

namespace stan {
namespace math {

/**
 * Forward mode specialization of calculating the inverse of the matrix.
 *
 * @tparam T type of elements in the matrix
 *
 * @param m specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not square.
 */
template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
inverse(const EigMat& m) {
  using T = typename value_type_t<EigMat>::Scalar;
  constexpr int R = EigMat::RowsAtCompileTime;
  constexpr int C = EigMat::ColsAtCompileTime;

  check_square("inverse", "m", m);
  if (m.size() == 0) {
    return {};
  }

  Eigen::Matrix<T, R, C> m_deriv(m.rows(), m.cols());
  Eigen::Matrix<T, R, C> m_inv(m.rows(), m.cols());

  const Eigen::Ref<const plain_type_t<EigMat>>& m_ref = m;
  for (size_type j = 0; j < m.cols(); j++) {
    for (size_type i = 0; i < m.rows(); i++) {
      m_inv.coeffRef(i, j) = m_ref.coeff(i, j).val_;
      m_deriv.coeffRef(i, j) = m_ref.coeff(i, j).d_;
    }
  }

  m_inv = inverse(m_inv);
  m_deriv = -multiply(multiply(m_inv, m_deriv), m_inv);

  return to_fvar(m_inv, m_deriv);
}

}  // namespace math
}  // namespace stan
#endif
