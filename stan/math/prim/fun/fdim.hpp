#ifndef STAN_MATH_PRIM_FUN_FDIM_HPP
#define STAN_MATH_PRIM_FUN_FDIM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * Return the positive difference of the specified values (C++11).
 *
 * The function is defined by
 *
 * <code>fdim(x, y) = (x > y) ? (x - y) : 0</code>.
 *
 * @param x First value.
 * @param y Second value.
 * @return max(x- y, 0)
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
inline double fdim(T1 x, T2 y) {
  using std::fdim;
  return fdim(x, y);
}

/**
 * Enables the vectorized application of the fdim function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Fdim function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto fdim(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return fdim(c, d); });
}

}  // namespace math
}  // namespace stan
#endif
