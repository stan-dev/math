#ifndef STAN_MATH_PRIM_FUN_MODIFIED_BESSEL_FIRST_KIND_HPP
#define STAN_MATH_PRIM_FUN_MODIFIED_BESSEL_FIRST_KIND_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{modified\_bessel\_first\_kind}(v, z) =
   \begin{cases}
     I_v(z) & \mbox{if } -\infty\leq z \leq \infty \\[6pt]
     \textrm{error} & \mbox{if } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{modified\_bessel\_first\_kind}(v, z)}{\partial z} =
   \begin{cases}
     \frac{\partial\, I_v(z)}{\partial z} & \mbox{if } -\infty\leq z\leq \infty
 \\[6pt] \textrm{error} & \mbox{if } z = \textrm{NaN} \end{cases} \f]

   \f[
     {I_v}(z) = \left(\frac{1}{2}z\right)^v\sum_{k=0}^\infty
 \frac{\left(\frac{1}{4}z^2\right)^k}{k!\Gamma(v+k+1)} \f]

   \f[
   \frac{\partial \, I_v(z)}{\partial z} = I_{v-1}(z)-\frac{v}{z}I_v(z)
   \f]
 *
 */
template <typename T2, require_arithmetic_t<T2>* = nullptr>
inline T2 modified_bessel_first_kind(int v, const T2 z) {
  check_not_nan("modified_bessel_first_kind", "z", z);

  return boost::math::cyl_bessel_i(v, z);
}

/**
 * This function exists because when z is of type integer,
 * cyl_bessel_i(v, z) returns an integer. This
 * results in overflow when the function value is large.
 */
inline double modified_bessel_first_kind(int v, int z) {
  check_not_nan("modified_bessel_first_kind", "z", z);

  return boost::math::cyl_bessel_i(v, static_cast<double>(z));
}

/**
 * Enables the vectorised application of the modified_bessel_first_kind
 * function, when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return modified_bessel_first_kind function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto modified_bessel_first_kind(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return modified_bessel_first_kind(c, d);
  });
}

}  // namespace math
}  // namespace stan

#endif
