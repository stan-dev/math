#ifndef STAN_MATH_PRIM_ARR_META_GET_HPP
#define STAN_MATH_PRIM_ARR_META_GET_HPP

#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <cstdlib>
#include <vector>

namespace stan {
/**
 * Returns the n-th element of the provided std::vector.
 *
 * @param x input vector
 * @param n index of the element to return
 * @return n-th element of the input vector
 */
template <typename T, require_std_vector<T>...>
inline auto&& get(T&& x, size_t n) {
  return std::forward<decltype(x[n])>(x[n]);
}

}  // namespace stan
#endif
