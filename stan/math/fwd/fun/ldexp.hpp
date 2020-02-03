#ifndef STAN_MATH_FWD_FUN_LDEXP_HPP
#define STAN_MATH_FWD_FUN_LDEXP_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/ldexp.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the product of a (the significand) times
 * 2 to power b (the exponent).
 *
 * @tparam T inner type of the fvar
 * @param[in] a the significand
 * @param[in] b an integer that is the exponent
 * @return product of a times 2 to the power b
 */
template <typename T>
inline fvar<T> ldexp(const fvar<T>& a, int b) {
  return fvar<T>(ldexp(a.val_, b), ldexp(a.d_, b));
}

}  // namespace math
}  // namespace stan

#endif
