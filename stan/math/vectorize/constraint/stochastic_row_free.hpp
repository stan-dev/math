
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_STOCHASTIC_ROW_FREE_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_STOCHASTIC_ROW_FREE_HPP
#include <stan/math/prim/constraint/stochastic_row_free.hpp>
namespace stan {
namespace math {

/**
 * Overload that untransforms each Eigen matrix in a standard vector.
 *
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and dynamic rows
 * @param[in] y vector of rowwise simplex matrices each of size (N, K)
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto stochastic_row_free(const T& y) {
  return apply_vector_unary<T>::apply(
      y, [](auto&& v) { return stochastic_row_free(v); });
}


} // namespace math
} // namespace stan
#endif 

