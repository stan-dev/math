
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_COV_MATRIX_FREE_LKJ_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_COV_MATRIX_FREE_LKJ_HPP
#include <stan/math/prim/constraint/cov_matrix_free_lkj.hpp>
namespace stan {
namespace math {

/**
 * Overload of `cov_matrix_free_lkj()` to untransform each matrix
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase`.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto cov_matrix_free_lkj(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [](auto&& v) { return cov_matrix_free_lkj(v); });
}


} // namespace math
} // namespace stan
#endif 

