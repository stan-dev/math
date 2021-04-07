#ifndef STAN_MATH_PRIM_FUN_BESSEL_FIRST_KIND_HPP
#define STAN_MATH_PRIM_FUN_BESSEL_FIRST_KIND_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{bessel\_first\_kind}(v, x) =
   \begin{cases}
     J_v(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{error} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{bessel\_first\_kind}(v, x)}{\partial x} =
   \begin{cases}
     \frac{\partial\, J_v(x)}{\partial x} & \mbox{if } -\infty\leq x\leq \infty
 \\[6pt] \textrm{error} & \mbox{if } x = \textrm{NaN} \end{cases} \f]

   \f[
   J_v(x)=\left(\frac{1}{2}x\right)^v
   \sum_{k=0}^\infty \frac{\left(-\frac{1}{4}x^2\right)^k}{k!\, \Gamma(v+k+1)}
   \f]

   \f[
   \frac{\partial \, J_v(x)}{\partial x} = \frac{v}{x}J_v(x)-J_{v+1}(x)
   \f]
 *
 */
template <typename T2, require_arithmetic_t<T2>* = nullptr>
inline T2 bessel_first_kind(int v, const T2 z) {
  check_not_nan("bessel_first_kind", "z", z);
  return boost::math::cyl_bessel_j(v, z);
}

/**
 * Enables the vectorised application of the bessel first kind function, when
 * the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Bessel first kind function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_not_var_matrix_t<T2>* = nullptr>
inline auto bessel_first_kind(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return bessel_first_kind(c, d);
  });
}

}  // namespace math
}  // namespace stan
#endif
