#ifndef STAN_MATH_PRIM_FUN_FMAX_HPP
#define STAN_MATH_PRIM_FUN_FMAX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the greater of the two specified arguments.  If one is
 * not-a-number, return the other.
 *
 * @param x First argument.
 * @param y Second argument.
 * @return maximum of x or y and if one is NaN return the other
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
inline double fmax(T1 x, T2 y) {
  using std::fmax;
  return fmax(x, y);
}

/**
 * Enables the vectorized application of the fmax function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return fmax function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto fmax(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return fmax(c, d); });
}

}  // namespace math
}  // namespace stan
#endif
