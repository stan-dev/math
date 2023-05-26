#ifndef STAN_MATH_PRIM_FUN_FMOD_HPP
#define STAN_MATH_PRIM_FUN_FMOD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the floating-point remainder of
 * numerator/denominator (rounded towards zero).
 *
 * @tparam T1 type of the numerator (must be arithmetic)
 * @tparam T2 type of the denominator (must be arithmetic)
 * @param a the numerator
 * @param b the denominator
 * @return the floating-point remainder of a/b.
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
inline double fmod(const T1& a, const T2& b) {
  return std::fmod(a, b);
}

/**
 * Enables the vectorized application of the fmod function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return fmod function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto fmod(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    using std::fmod;
    return fmod(c, d);
  });
}

}  // namespace math
}  // namespace stan
#endif
