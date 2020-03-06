#ifndef STAN_MATH_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_PRIM_FUN_INV_SQUARE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/inv_square.hpp>

namespace stan {
namespace math {

/**
 * Returns 1 / square(x).
 *
 * @tparam T type of container
 * @param x container
 * @return 1 / square(x) of each value in x.
 */
template <typename T>
inline auto inv_square(const T& x) {
  return inv(square(x));
}

}  // namespace math
}  // namespace stan

#endif
