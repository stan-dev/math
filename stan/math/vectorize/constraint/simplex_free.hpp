
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_SIMPLEX_FREE_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_SIMPLEX_FREE_HPP
#include <stan/math/prim/constraint/simplex_free.hpp>
namespace stan {
namespace math {

/**
 * Overload of `simplex_free()` to untransform each Eigen vector
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase` with compile time rows or columns equal to 1.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto simplex_free(const T& x) {
  return apply_vector_unary<T>::apply(x,
                                      [](auto&& v) { return simplex_free(v); });
}


} // namespace math
} // namespace stan
#endif 

