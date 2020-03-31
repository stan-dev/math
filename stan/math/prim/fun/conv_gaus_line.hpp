#ifndef STAN_MATH_PRIM_FUN_CONV_GAUS_LINE
#define STAN_MATH_PRIM_FUN_CONV_GAUS_LINE

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/normal_cdf.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/*
evaluate the derivative of conv_gaus_line with respect to x
*/
double der_conv_gaus_line(double t0, double t1, double a, double b, double x0, 
			  double sig2) {
  using stan::math::normal_cdf;
  double pi = stan::math::pi();
  using std::exp;
  using std::pow;
  using std::sqrt;
  double sig = sqrt(sig2);
  double y;

  double alpha = sqrt(2 * pi * sig2);
  y = (a * x0 + b) / alpha * (-exp(-pow(t1 - x0, 2) / (2 * sig2))
			      + exp(-pow(t0 - x0, 2) / (2 * sig2)));
  y +=  a * (normal_cdf(t1, x0, sig) - normal_cdf(t0, x0, sig));
  y -= a / alpha * ((t1-x0)*exp(-pow(t1 - x0, 2) / (2 * sig2))
		    -(t0-x0)*exp(-pow(t0 - x0, 2) / (2 * sig2)));
  return y;
}


/*
evaluate the integral

                    | t1
(2*pi*sig2)**-(1/2) |   (a*t + b) * exp(-(t - x0)^2 / (2*sig2)) dt
                    | t0

*/
template <typename Tx>
double conv_gaus_line(double t0, double t1, double a, double b, Tx const& x, 
		      double sig2) {
  using stan::math::normal_cdf;
  double pi = stan::math::pi();
  using std::exp;
  using std::pow;
  using std::sqrt;
  double sig = sqrt(sig2);

  double y = (a * value_of(x) + b) * (normal_cdf(t1, value_of(x), sig) 
			     - normal_cdf(t0, value_of(x), sig));
  y += - a * sig2 / sqrt(2 * pi * sig2) * 
    (exp(-pow(t1 - value_of(x), 2)/ (2 * sig2)) 
     - exp(-pow(t0 - value_of(x), 2) / (2 * sig2)));

  return y;
}

}  // namespace math
}  // namespace stan

#endif
