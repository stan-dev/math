#ifndef STAN_MATH_PRIM_FUN_COV_MATRIX_FREE_HPP
#define STAN_MATH_PRIM_FUN_COV_MATRIX_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The covariance matrix derived from the symmetric view of the
 * lower-triangular view of the K by K specified matrix is freed
 * to return a vector of size K + (K choose 2).
 *
 * This is the inverse of the <code>cov_matrix_constrain()</code>
 * function so that for any finite vector <code>x</code> of size K
 * + (K choose 2),
 *
 * <code>x == cov_matrix_free(cov_matrix_constrain(x, K))</code>.
 *
 * In order for this round-trip to work (and really for this
 * function to work), the symmetric view of its lower-triangular
 * view must be positive definite.
 *
 * @tparam T type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param y Matrix of dimensions K by K such that he symmetric
 * view of the lower-triangular view is positive definite.
 * @return Vector of size K plus (K choose 2) in (-inf, inf)
 * that produces
 * @throw std::domain_error if <code>y</code> is not square,
 * has zero dimensionality, or has a non-positive diagonal element.
 */
template <typename T, require_eigen_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1> cov_matrix_free(const T& y) {
  const auto& y_ref = to_ref(y);
  check_square("cov_matrix_free", "y", y_ref);
  check_nonzero_size("cov_matrix_free", "y", y_ref);

  using std::log;
  int K = y_ref.rows();
  check_positive("cov_matrix_free", "y", y_ref.diagonal());
  Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1> x((K * (K + 1)) / 2);
  // FIXME: see Eigen LDLT for rank-revealing version -- use that
  // even if less efficient?
  Eigen::LLT<plain_type_t<T>> llt(y_ref.rows());
  llt.compute(y_ref);
  plain_type_t<T> L = llt.matrixL();
  int i = 0;
  for (int m = 0; m < K; ++m) {
    x.segment(i, m) = L.row(m).head(m);
    i += m;
    x.coeffRef(i++) = log(L.coeff(m, m));
  }
  return x;
}

/**
 * Overload of `cov_matrix_free()` to untransform each matrix
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase`.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto cov_matrix_free(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](auto&& v) { return cov_matrix_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
