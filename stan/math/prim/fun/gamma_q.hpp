#ifndef STAN_MATH_PRIM_FUN_GAMMA_Q_HPP
#define STAN_MATH_PRIM_FUN_GAMMA_Q_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace stan {
namespace math {

/**
 *
   \f[
   \mbox{gamma\_q}(a, z) =
   \begin{cases}
     \textrm{error} & \mbox{if } a\leq 0 \textrm{ or } z < 0\\
     Q(a, z) & \mbox{if } a > 0, z \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } a = \textrm{NaN or } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{gamma\_q}(a, z)}{\partial a} =
   \begin{cases}
     \textrm{error} & \mbox{if } a\leq 0 \textrm{ or } z < 0\\
     \frac{\partial\, Q(a, z)}{\partial a} & \mbox{if } a > 0, z \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } a = \textrm{NaN or } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{gamma\_q}(a, z)}{\partial z} =
   \begin{cases}
     \textrm{error} & \mbox{if } a\leq 0 \textrm{ or } z < 0\\
     \frac{\partial\, Q(a, z)}{\partial z} & \mbox{if } a > 0, z \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } a = \textrm{NaN or } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   Q(a, z)=\frac{1}{\Gamma(a)}\int_z^\infty t^{a-1}e^{-t}dt
   \f]

   \f[
   \frac{\partial \, Q(a, z)}{\partial a} =
 -\frac{\Psi(a)}{\Gamma^2(a)}\int_z^\infty t^{a-1}e^{-t}dt
   + \frac{1}{\Gamma(a)}\int_z^\infty (a-1)t^{a-2}e^{-t}dt
   \f]

   \f[
   \frac{\partial \, Q(a, z)}{\partial z} = -\frac{z^{a-1}e^{-z}}{\Gamma(a)}
   \f]
   * @throws domain_error if x is at pole
 */
inline double gamma_q(double x, double a) { return boost::math::gamma_q(x, a); }

/**
 * Enables the vectorized application of the gamma_q function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return gamma_q function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto gamma_q(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return gamma_q(c, d); });
}

}  // namespace math
}  // namespace stan
#endif
