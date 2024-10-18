
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_SUM_TO_ZERO_FREE_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_SUM_TO_ZERO_FREE_HPP
#include <stan/math/prim/constraint/sum_to_zero_free.hpp>
namespace stan {
namespace math {

/**
 * Overload of `sum_to_zero_free()` to untransform each Eigen vector
 * in a standard vector.
 * @tparam T A standard vector with with a `value_type` which inherits from
 *  `Eigen::MatrixBase` with compile time rows or columns equal to 1.
 * @param z The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto sum_to_zero_free(const T& z) {
  return apply_vector_unary<T>::apply(
      z, [](auto&& v) { return sum_to_zero_free(v); });
}


} // namespace math
} // namespace stan
#endif 

