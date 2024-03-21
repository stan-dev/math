#ifndef STAN_MATH_PRIM_FUN_GAMMA_P_HPP
#define STAN_MATH_PRIM_FUN_GAMMA_P_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the normalized, lower-incomplete gamma function
 * applied to the specified argument.
 *
 * <p>This function is defined, including error conditions, as follows
   \f[
   \mbox{gamma\_p}(a, z) =
   \begin{cases}
     \textrm{error} & \mbox{if } a\leq 0 \textrm{ or } z < 0\\
     P(a, z) & \mbox{if } a > 0, z \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } a = \textrm{NaN or } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{gamma\_p}(a, z)}{\partial a} =
   \begin{cases}
     \textrm{error} & \mbox{if } a\leq 0 \textrm{ or } z < 0\\
     \frac{\partial\, P(a, z)}{\partial a} & \mbox{if } a > 0, z \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } a = \textrm{NaN or } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{gamma\_p}(a, z)}{\partial z} =
   \begin{cases}
     \textrm{error} & \mbox{if } a\leq 0 \textrm{ or } z < 0\\
     \frac{\partial\, P(a, z)}{\partial z} & \mbox{if } a > 0, z \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } a = \textrm{NaN or } z = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   P(a, z)=\frac{1}{\Gamma(a)}\int_0^zt^{a-1}e^{-t}dt
   \f]

   \f[
   \frac{\partial \, P(a, z)}{\partial a} =
 -\frac{\Psi(a)}{\Gamma^2(a)}\int_0^zt^{a-1}e^{-t}dt
   + \frac{1}{\Gamma(a)}\int_0^z (a-1)t^{a-2}e^{-t}dt
   \f]

   \f[
   \frac{\partial \, P(a, z)}{\partial z} = \frac{z^{a-1}e^{-z}}{\Gamma(a)}
   \f]
   *
   * @param z first argument
   * @param a second argument
   * @return value of the normalized, lower-incomplete gamma function
   * applied to z and a
   * @throws std::domain_error if either argument is not positive or
   * if z is at a pole of the function
 */
inline double gamma_p(double z, double a) {
  if (is_nan(z)) {
    return not_a_number();
  }
  if (is_nan(a)) {
    return not_a_number();
  }
  check_positive("gamma_p", "first argument (z)", z);
  check_nonnegative("gamma_p", "second argument (a)", a);
  return boost::math::gamma_p(z, a, boost_policy_t<>());
}

/**
 * Enables the vectorized application of the gamma_p function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return gamma_p function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto gamma_p(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return gamma_p(c, d); });
}

}  // namespace math
}  // namespace stan
#endif
