#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_integrate_1d_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_integrate_1d_HPP

#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/functor/de_integrator.hpp>
#include <functional>
#include <cmath>
#include <ostream>
#include <vector>

namespace stan {

namespace math {

/**
 * Return the numerical integral of a function f.
 *
 * The numerical integration algorithm used is the Tanh-sinh
 * quadrature method (also known as Double Exponential
 * Transformation) which was proposed by Hidetosi Takahasi and
 * Masatake Mori in 1974.
 *
 * The implementation of integration used is given John D. Cook. See
 * www.codeproject.com/kb/recipes/fastnumericalintegration.aspx
 * www.johndcook.com/blog/double_exponential_integration/
 *
 * The signature for the function to be integrated is:
 * double (double x, std::vector<double> params, std::vector<double> x_r,
 *   std::vector<int> x_i, std::ostream& msgs)
 *
 * It should return the value of the function evaluated at x. Any errors should
 * be printed to the msgs stream.
 *
 * Such implementation assumes f to be smooth with no discontinuity
 * in the function nor in any of its derivatives.
 *
 * The integration is terminated when the absolute value of the current estimate
 * of the integral minus the previous estimate of the integral is less than
 * tolerance.
 *
 * @tparam T Type of f
 * @param f a functor to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param tolerance target absolute error
 * @param param additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param msgs stream
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate_1d(const F& f, const double a, const double b,
                           const std::vector<double>& param,
                           const std::vector<double>& x_r,
                           const std::vector<int>& x_i, std::ostream& msgs,
                           const double tolerance = 1e-6) {
  using std::placeholders::_1;
  static const char* function = "integrate_1d";

  stan::math::check_finite(function, "lower limit", a);
  stan::math::check_finite(function, "upper limit", b);
  stan::math::check_less_or_equal(function, "lower limit", a, b);

  double val_
      = de_integrator(std::bind<double>(f, _1, param, x_r, x_i, std::ref(msgs)),
                      a, b, tolerance);

  return val_;
}

}  // namespace math

}  // namespace stan

#endif
