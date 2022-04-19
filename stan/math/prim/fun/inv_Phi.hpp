#ifndef STAN_MATH_PRIM_FUN_INV_PHI_HPP
#define STAN_MATH_PRIM_FUN_INV_PHI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
/**
 *  The largest integer that protects against floating point errors
 * for the inv_Phi function. The value was found by finding the largest
 * integer that passed the unit tests for accuracy when the input into inv_Phi
 * is near 1.
 */
const int BIGINT = 2000000000;

/**
 * The inverse of the unit normal cumulative distribution function.
 *
 * @param p argument between 0 and 1 inclusive
 * @return Real value of the inverse cdf for the standard normal distribution.
 */
inline double inv_Phi_lambda(double p) {
  check_bounded("inv_Phi", "Probability variable", p, 0, 1);

  if (p < 8e-311) {
    return NEGATIVE_INFTY;
  }
  if (p == 1) {
    return INFTY;
  }

  static const double a[8]
      = {3.3871328727963666080e+00, 1.3314166789178437745e+02,
         1.9715909503065514427e+03, 1.3731693765509461125e+04,
         4.5921953931549871457e+04, 6.7265770927008700853e+04,
         3.3430575583588128105e+04, 2.5090809287301226727e+03};
  static const double b[7]
      = {4.2313330701600911252e+01, 6.8718700749205790830e+02,
         5.3941960214247511077e+03, 2.1213794301586595867e+04,
         3.9307895800092710610e+04, 2.8729085735721942674e+04,
         5.2264952788528545610e+03};
  static const double c[8]
      = {1.42343711074968357734e+00, 4.63033784615654529590e+00,
         5.76949722146069140550e+00, 3.64784832476320460504e+00,
         1.27045825245236838258e+00, 2.41780725177450611770e-01,
         2.27238449892691845833e-02, 7.74545014278341407640e-04};
  static const double d[7]
      = {2.05319162663775882187e+00, 1.67638483018380384940e+00,
         6.89767334985100004550e-01, 1.48103976427480074590e-01,
         1.51986665636164571966e-02, 5.47593808499534494600e-04,
         1.05075007164441684324e-09};
  static const double e[8]
      = {6.65790464350110377720e+00, 5.46378491116411436990e+00,
         1.78482653991729133580e+00, 2.96560571828504891230e-01,
         2.65321895265761230930e-02, 1.24266094738807843860e-03,
         2.71155556874348757815e-05, 2.01033439929228813265e-07};
  static const double f[7]
      = {5.99832206555887937690e-01, 1.36929880922735805310e-01,
         1.48753612908506148525e-02, 7.86869131145613259100e-04,
         1.84631831751005468180e-05, 1.42151175831644588870e-07,
         2.04426310338993978564e-15};

  double q = p - 0.5;
  double r;
  double val;

  if (std::fabs(q) <= .425) {
    r = .180625 - square(q);
    return q
           * (((((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3]) * r
                + a[2])
                   * r
               + a[1])
                  * r
              + a[0])
           / (((((((b[6] * r + b[5]) * r + b[4]) * r + b[3]) * r + b[2]) * r
                + b[1])
                   * r
               + b[0])
                  * r
              + 1.0);
  } else {
    r = q < 0 ? p : 1 - p;

    if (r <= 0)
      return 0;

    r = std::sqrt(-std::log(r));

    if (r <= 5.0) {
      r += -1.6;
      val = (((((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3]) * r
               + c[2])
                  * r
              + c[1])
                 * r
             + c[0])
            / (((((((d[6] * r + d[5]) * r + d[4]) * r + d[3]) * r + d[2]) * r
                 + d[1])
                    * r
                + d[0])
                   * r
               + 1.0);
    } else {
      r -= 5.0;
      val = (((((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3]) * r
               + e[2])
                  * r
              + e[1])
                 * r
             + e[0])
            / (((((((f[6] * r + f[5]) * r + f[4]) * r + f[3]) * r + f[2]) * r
                 + f[1])
                    * r
                + f[0])
                   * r
               + 1.0);
    }
    if (q < 0.0)
      return -val;
  }
  return val;
}
}  // namespace internal

/**
 * Return the value of the inverse standard normal cumulative distribution
 * function at the specified argument.
 *
 * The precision is at or better than 1.5e-15 for values between 0.0000001 he
 * largest integer that protects against floating point errors for the inv_Phi
 * function. The value was found by finding the largest integer that passed the
 * unit tests for accuracy when the input into inv_Phi is near 1.
 *
 * @param p argument between 0 and 1 inclusive
 * @return real value of the inverse cdf for the standard normal distribution
 */
inline double inv_Phi(double p) {
  return p >= 0.9999 ? -internal::inv_Phi_lambda(
             (internal::BIGINT - internal::BIGINT * p) / internal::BIGINT)
                     : internal::inv_Phi_lambda(p);
}

/**
 * Structure to wrap inv_Phi() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable in range [0, 1]
 * @return Inverse unit normal CDF of x.
 * @throw std::domain_error if x is not between 0 and 1.
 */
struct inv_Phi_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_Phi(x);
  }
};

/**
 * Vectorized version of inv_Phi().
 *
 * @tparam T type of container
 * @param x variables in range [0, 1]
 * @return Inverse unit normal CDF of each value in x.
 * @throw std::domain_error if any value is not between 0 and 1.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto inv_Phi(const T& x) {
  return apply_scalar_unary<inv_Phi_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
