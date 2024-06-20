#ifndef STAN_MATH_PRIM_FUN_MAX_SIZE_HPP
#define STAN_MATH_PRIM_FUN_MAX_SIZE_HPP

#include <stan/math/prim/fun/size.hpp>
#include <cstdint>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Calculate the size of the largest input.
 * @tparam T1 type of the first input
 * @tparam Ts types of the other inputs
 * @param x1 first input
 * @param xs other inputs
 * @return the size of the largest input
 */
template <typename T1, typename... Ts>
inline int64_t max_size(const T1& x1, const Ts&... xs) {
  return std::max({stan::math::size(x1), stan::math::size(xs)...});
}

}  // namespace math
}  // namespace stan
#endif
