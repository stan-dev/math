#ifndef STAN_MATH_PRIM_FUN_MODIFIED_BESSEL_SECOND_KIND_HPP
#define STAN_MATH_PRIM_FUN_MODIFIED_BESSEL_SECOND_KIND_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{modified\_bessel\_second\_kind}(v, z) =
   \begin{cases}
     \textrm{error} & \mbox{if } z \leq 0 \\
     K_v(z) & \mbox{if } z > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{modified\_bessel\_second\_kind}(v, z)}{\partial z} =
   \begin{cases}
     \textrm{error} & \mbox{if } z \leq 0 \\
     \frac{\partial\, K_v(z)}{\partial z} & \mbox{if } z > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   {K_v}(z)
   =
   \frac{\pi}{2}\cdot\frac{I_{-v}(z) - I_{v}(z)}{\sin(v\pi)}
   \f]

   \f[
   \frac{\partial \, K_v(z)}{\partial z} = -\frac{v}{z}K_v(z)-K_{v-1}(z)
   \f]
 *
 */
template <typename T2, require_arithmetic_t<T2>* = nullptr>
inline T2 modified_bessel_second_kind(int v, const T2 z) {
  return boost::math::cyl_bessel_k(v, z);
}

/**
 * Enables the vectorized application of the modified_bessel_second_kind
 * function, when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return modified_bessel_second_kind function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto modified_bessel_second_kind(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return modified_bessel_second_kind(c, d);
  });
}

}  // namespace math
}  // namespace stan

#endif
