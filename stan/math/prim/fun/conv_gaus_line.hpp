#ifndef STAN_MATH_PRIM_FUN_CONV_GAUS_LINE
#define STAN_MATH_PRIM_FUN_CONV_GAUS_LINE

#include <cmath>
#include <vector>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/normal_cdf.hpp>

namespace stan {
namespace math {


/** \ingroup prob_dists
 * Evaluates the derivative of the convolution of a line with a Gaussian
 * kernel on an interval.
 *
 * \f$\frac{\partial}{\partial x} \int_{t_0}^{t_1} (at + b) \f$
 * \f$ e^{\frac{-(t-x)^2}{2\sigma^2}} dt \f$
 *
 * @param t_0 lower integration bound
 * @param t_1 upper integration bound
 * @param a coefficient of t in line
 * @param b constant in line
 * @param x0 point at which convolution is evaluated 
 * @param sig2 variance of the Gaussian kernel
 * @return The value of the derivative
 */
double der_conv_gaus_line(double t0, double t1, double a, double b, double x0,
                          double sig2) {
  using stan::math::normal_cdf;
  using std::exp;
  using std::pow;
  using std::sqrt;
  const double pi = stan::math::pi();
  const double sig = sqrt(sig2);
  const double alpha = sqrt(2 * pi * sig2);

  double y = (a * x0 + b) / alpha
      * (-exp(-pow(t1 - x0, 2) / (2 * sig2))
         + exp(-pow(t0 - x0, 2) / (2 * sig2)));
  y += a * (normal_cdf(t1, x0, sig) - normal_cdf(t0, x0, sig));
  y -= a / alpha
       * ((t1 - x0) * exp(-pow(t1 - x0, 2) / (2 * sig2))
          - (t0 - x0) * exp(-pow(t0 - x0, 2) / (2 * sig2)));
  return y;
}


/** \ingroup prob_dists
 * Evaluate the convolution of a line with a Gaussian kernel on an interval.
 *
 * \f$\int_{t_0}^{t_1} (at + b) e^{\frac{-(t-x)^2}{2\sigma^2}} dt \f$
 *
 * @param t_0 lower integration bound
 * @param t_1 upper integration bound
 * @param a coefficient of t in line
 * @param b constant in line
 * @param x0 point at which convolution is evaluated 
 * @param sig2 variance of the Gaussian kernel
 * @return The value of the derivative
 */
template <typename Tx>
double conv_gaus_line(double t0, double t1, double a, double b, Tx const& x,
                      double sig2) {
  using stan::math::normal_cdf;
  using std::exp;
  using std::pow;
  using std::sqrt;
  const double pi = stan::math::pi();
  const double sig = sqrt(sig2);

  double y
      = (a * value_of(x) + b)
        * (normal_cdf(t1, value_of(x), sig) - normal_cdf(t0, value_of(x), sig));
  y += -a * sig2 / sqrt(2 * pi * sig2)
       * (exp(-pow(t1 - value_of(x), 2) / (2 * sig2))
          - exp(-pow(t0 - value_of(x), 2) / (2 * sig2)));
  return y;
}

}  // namespace math
}  // namespace stan

#endif
