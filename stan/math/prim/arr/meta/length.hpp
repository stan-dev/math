#ifndef STAN_MATH_PRIM_ARR_META_LENGTH_HPP
#define STAN_MATH_PRIM_ARR_META_LENGTH_HPP

#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <cstdlib>
#include <vector>

namespace stan {
/**
 * Returns the length of the provided std::vector.
 *
 * @param x input vector
 * @tparam T type of the elements in the vector
 * @return the length of the input vector
 */
template <typename T, require_std_vector<T>...>
auto&& length(T&& x) {
  return std::forward<decltype(x.size())>(x.size());
}
}  // namespace stan
#endif
