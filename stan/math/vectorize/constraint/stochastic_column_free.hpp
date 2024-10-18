
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_STOCHASTIC_COLUMN_FREE_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_STOCHASTIC_COLUMN_FREE_HPP
#include <stan/math/prim/constraint/stochastic_column_free.hpp>
namespace stan {
namespace math {

/**
 * Overload that untransforms each Eigen matrix in a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic rows
 * @param[in] y vector of columnwise stochastic matrix of size (N, K)
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_column_free(const T& y) {
  return apply_vector_unary<T>::apply(
      y, [](auto&& v) { return stochastic_column_free(v); });
}


} // namespace math
} // namespace stan
#endif 

