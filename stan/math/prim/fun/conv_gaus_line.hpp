#ifndef STAN_MATH_PRIM_FUN_CONV_GAUS_LINE
#define STAN_MATH_PRIM_FUN_CONV_GAUS_LINE

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/normal_cdf.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Evaluate the convolution of a line with a Gaussian kernel on an interval.
 *
 * \f$\int_{t_0}^{t_1} (at + b) e^{\frac{-(t-x)^2}{2\sigma^2}} dt \f$
 *
 * @param t0 lower integration bound
 * @param t1 upper integration bound
 * @param a coefficient of t in line
 * @param b constant in line
 * @param x point at which convolution is evaluated
 * @param sig2 variance of the Gaussian kernel
 * @return The value of the derivative
 */
template <typename Tx>
inline double conv_gaus_line(double t0, double t1, double a, double b,
                             const Tx& x, double sig2) {
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


/**
 * Evaluate the convolution of a piecewise linear function with a 
 * Gaussian kernel.
 *
 * \f$\int_{x_0}^{x_1} (a_1t + b_1) e^{\frac{-(t-x)^2}{2\sigma^2}} dt \f$
 * \f$ + ... + $\int_{x_{n-1}}^{x_{n}} (a_{n-1}t + b_{n-1}) \f$
 * \f$ e^{\frac{-(t-x)^2}{2\sigma^2}} dt \f$
 *
 * @param x point at which convolution is evaluated
 * @param xs increasing vector with endpoints of each piecewise linear function
 * @param params vector containing slopes and intercepts of peicewise linear function and width of Gaussian kernel
 * @param ind_start
 * @param ind_end
 * @return The value of the derivative
 */
template <typename Tx>
inline double conv_gaus_line_sum(const Tx& x, std::vector<double> xs,
				 std::vector<double> params, 
				 int ind_start, 
				 int ind_end) {
  using stan::math::normal_cdf;
  using std::exp;
  using std::pow;
  using std::sqrt;
  const double pi = stan::math::pi();
  int n = xs.size();
  double sig2 = params[2*n-2];
  double sig = std::sqrt(sig2);

  double y=0;
  for (int i = ind_start; i < ind_end; i++) {
    y += (params[i] * value_of(x) + params[n-1 + i])
      * (normal_cdf(xs[i+1], value_of(x), sig) - normal_cdf(xs[i], value_of(x), sig));
    y += -params[i] * sig2 / sqrt(2 * pi * sig2)
      * (exp(-pow(xs[i+1] - value_of(x), 2) / (2 * sig2))
	 - exp(-pow(xs[i] - value_of(x), 2) / (2 * sig2)));
  }
  return y;
}

}  // namespace math
}  // namespace stan

#endif
