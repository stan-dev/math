#ifndef STAN_MATH_PRIM_FUN_ATAN2_HPP
#define STAN_MATH_PRIM_FUN_ATAN2_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Computes the arc tangent of y/x using the signs of
 * arguments to determine the correct quadrant.
 *
 * @tparam T1 type of first argument (must be arithmetic)
 * @tparam T2 type of second argument (must be arithmetic)
 * @param y first argument
 * @param x second argument
 * @return arctangent of `y / x`
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
double atan2(T1 y, T2 x) {
  return std::atan2(y, x);
}

/**
 * Enables the vectorised application of the atan2 function, when
 * the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Returns the atan2 function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_all_not_var_matrix_t<T1, T2>* = nullptr>
inline auto atan2(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [](const auto& c, const auto& d) { return atan2(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
