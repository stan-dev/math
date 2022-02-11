#ifndef STAN_MATH_PRIM_FUN_BESSEL_SECOND_KIND_HPP
#define STAN_MATH_PRIM_FUN_BESSEL_SECOND_KIND_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{bessel\_second\_kind}(v, x) =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0 \\
     Y_v(x) & \mbox{if } x > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{bessel\_second\_kind}(v, x)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0 \\
     \frac{\partial\, Y_v(x)}{\partial x} & \mbox{if } x > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   Y_v(x)=\frac{J_v(x)\cos(v\pi)-J_{-v}(x)}{\sin(v\pi)}
   \f]

   \f[
   \frac{\partial \, Y_v(x)}{\partial x} = \frac{v}{x}Y_v(x)-Y_{v+1}(x)
   \f]
 *
 */
template <typename T2, require_arithmetic_t<T2>* = nullptr>
inline T2 bessel_second_kind(int v, const T2 z) {
  return boost::math::cyl_neumann(v, z);
}

/**
 * Enables the vectorised application of the bessel second kind function, when
 * the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Bessel second kind function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto bessel_second_kind(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return bessel_second_kind(c, d);
  });
}

}  // namespace math
}  // namespace stan

#endif
