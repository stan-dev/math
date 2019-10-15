#ifndef STAN_MATH_PRIM_SCAL_FUN_GAMMA_P_HPP
#define STAN_MATH_PRIM_SCAL_FUN_GAMMA_P_HPP

#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>

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
  return boost::math::gamma_p(z, a, boost_policy_t());
}

}  // namespace math
}  // namespace stan
#endif
