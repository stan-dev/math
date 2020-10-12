#ifndef STAN_MATH_PRIM_ERR_IS_POS_DEFINITE_HPP
#define STAN_MATH_PRIM_ERR_IS_POS_DEFINITE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/is_positive.hpp>
#include <stan/math/prim/err/is_not_nan.hpp>
#include <stan/math/prim/err/is_symmetric.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the matrix is square or if the matrix has
 * non-zero size, or if the matrix is symmetric, or if it is positive
 * definite, or if no element is <code>NaN</code>.
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square or if the matrix has non-0
 *   size, or if the matrix is symmetric, or if it is positive definite, or if
 *   no element is <code>NaN</code>
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline bool is_pos_definite(const EigMat& y) {
  const auto& y_ref = to_ref(y);
  if (!is_symmetric(y_ref)) {
    return false;
  }
  if (!is_positive(y_ref.rows())) {
    return false;
  }
  if (y_ref.rows() == 1 && !(y_ref(0, 0) > CONSTRAINT_TOLERANCE)) {
    return false;
  }
  Eigen::LDLT<Eigen::MatrixXd> cholesky = value_of_rec(y_ref).ldlt();
  if (cholesky.info() != Eigen::Success || !cholesky.isPositive()
      || (cholesky.vectorD().array() <= 0.0).any()) {
    return false;
  }
  return is_not_nan(y_ref);
}

/**
 * Return <code>true</code> if the matrix is positive definite.
 * @tparam Derived Derived type of the Eigen::LDLT transform, requires
 *   class method <code>.info()</code> and <code>.isPositive()</code>
 *   and <code>.vectorD().array()</code>
 * @param cholesky Eigen::LDLT to test, whose progenitor
 * must not have any NaN elements
 * @return <code>true</code> if the matrix is positive definite
 */
template <typename Derived>
inline bool is_pos_definite(const Eigen::LDLT<Derived>& cholesky) {
  return cholesky.info() == Eigen::Success && cholesky.isPositive()
         && (cholesky.vectorD().array() > 0.0).all();
}

/**
 * Return <code>true</code> if diagonal of the L matrix is positive.
 * @tparam Derived Derived type of the Eigen::LLT transform, requires
 *   class method <code>.info()</code> and
 *   <code>.matrixLLT().diagonal().array()</code>
 * @param cholesky Eigen::LLT to test, whose progenitor
 * must not have any NaN elements
 * @return <code>true</code> if diagonal of the L matrix is positive
 */
template <typename Derived>
inline bool is_pos_definite(const Eigen::LLT<Derived>& cholesky) {
  return cholesky.info() == Eigen::Success
         && (cholesky.matrixLLT().diagonal().array() > 0.0).all();
}

}  // namespace math
}  // namespace stan
#endif
