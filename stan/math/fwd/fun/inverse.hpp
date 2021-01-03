#ifndef STAN_MATH_FWD_FUN_INVERSE_HPP
#define STAN_MATH_FWD_FUN_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/inverse.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
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
inline auto inverse(const EigMat& m) {
  constexpr int R = EigMat::RowsAtCompileTime;
  constexpr int C = EigMat::ColsAtCompileTime;

  check_square("inverse", "m", m);
  if (m.size() == 0) {
    return Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                         EigMat::ColsAtCompileTime>();
  }
  const auto& m_ref = to_ref(m);
  auto m_inv = inverse(m_ref.val()).eval();
  auto m_deriv = (-multiply(multiply(m_inv, m_ref.d().eval()), m_inv)).eval();
  Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                       EigMat::ColsAtCompileTime> ret(m_inv);
  ret.d() = m_deriv;
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
