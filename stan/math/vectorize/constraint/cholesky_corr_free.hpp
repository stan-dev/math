
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_CHOLESKY_CORR_FREE_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_CHOLESKY_CORR_FREE_HPP
#include <stan/math/prim/constraint/cholesky_corr_free.hpp>
namespace stan {
namespace math {

/**
 * Overload of `cholesky_corr_free()` to untransform each matrix
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase`.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto cholesky_corr_free(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](auto&& v) { return cholesky_corr_free(v); });
}


} // namespace math
} // namespace stan
#endif 

