#ifndef STAN_MATH_REV_FUN_CONV_GAUS_LINE
#define STAN_MATH_REV_FUN_CONV_GAUS_LINE

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/conv_gaus_line.hpp>

namespace stan {
namespace math {

namespace internal {

/** \ingroup prob_dists
 * Evaluates the derivative of the convolution of a line with a Gaussian
 * kernel on an interval.
 *
 * \f$\frac{\partial}{\partial x} \int_{t_0}^{t_1} (at + b) \f$
 * \f$ e^{\frac{-(t-x)^2}{2\sigma^2}} dt \f$
 *
 * @param t0 lower integration bound
 * @param t1 upper integration bound
 * @param a coefficient of t in line
 * @param b constant in line
 * @param x0 point at which convolution is evaluated
 * @param sig2 variance of the Gaussian kernel
 * @return The value of the derivative
 */
inline double der_conv_gaus_line(double t0, double t1, double a, double b,
                                 double x0, double sig2) {
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

class conv_gaus_line_vari : public op_v_vari {
  double t0_;
  double t1_;
  double a_;
  double b_;
  double sig2_;

 public:
  explicit conv_gaus_line_vari(double t0, double t1, double a, double b,
                               vari* avi, double sig2)
      : op_v_vari(conv_gaus_line(t0, t1, a, b, avi->val_, sig2), avi),
        t0_(t0),
        t1_(t1),
        a_(a),
        b_(b),
        sig2_(sig2) {}

  void chain() {
    avi_->adj_
        += adj_ * der_conv_gaus_line(t0_, t1_, a_, b_, avi_->val_, sig2_);
  }
};

}  // namespace internal

inline var conv_gaus_line(double t0, double t1, double a, double b,
                          const var& x, double sig2) {
  return var(new internal::conv_gaus_line_vari(t0, t1, a, b, x.vi_, sig2));
}

}  // namespace math
}  // namespace stan

#endif
