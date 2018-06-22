#ifndef STAN_MATH_REV_ARR_FUNCTOR_integrate_1d_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_integrate_1d_HPP

#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/arr/functor/integrate_1d.hpp>
#include <functional>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Calculate first derivative of f(x, param, std::ostream&)
 * with respect to the nth parameter. Uses nested reverse mode autodiff
 */
template <typename F>
inline double gradient_of_f(const F& f, const double x,
                            const std::vector<double>& theta_vals,
                            const std::vector<double>& x_r,
                            const std::vector<int>& x_i, size_t n,
                            std::ostream& msgs) {
  double gradient = 0.0;
  start_nested();
  std::vector<var> theta_var(theta_vals.size());
  try {
    for (size_t i = 0; i < theta_vals.size(); i++)
      theta_var[i] = theta_vals[i];
    f(x, theta_var, x_r, x_i, msgs).grad();
    gradient = theta_var[n].adj();
  } catch (const std::exception& e) {
    recover_memory_nested();
    throw;
  }
  recover_memory_nested();

  return gradient;
}

/**
 * Compute the integral of the single variable function f from a to b to within
 * a specified tolerance. a and b can be finite or infinite. a and b cannot be
 * parameters.
 *
 * f should be compatible with reverse mode autodiff and have the signature:
 *   var foo(double x, std::vector<var> theta, std::vector<double> x_r,
 *     std::vector<int> x_i, std::ostream& msgs)
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream.
 *
 * Boost decides the integration is converged when subsequent estimates of the
 * integral are less than tolerance * abs(integral). This means the tolerance is
 * relative to the actual scale of the integral. Integrals that cross zero are
 * split into two. In this case, each integral is separately integrated to the
 * given tolerance.
 *
 * @tparam T Type of f
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param theta additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param tolerance integrator tolerance passed to Boost quadrature
 * @return numeric integral of function f
 */
template <typename F>
inline var integrate_1d(const F& f, const double a, const double b,
                        const std::vector<var>& theta,
                        const std::vector<double>& x_r,
                        const std::vector<int>& x_i, std::ostream& msgs,
                        const double tolerance = 1e-6) {
  static const char* function = "integrate_1d";
  check_less_or_equal(function, "lower limit", a, b);

  if (a == b) {
    if (std::isinf(a))
      domain_error(function, "Integration endpoints are both", a, "", "");
    return var(0.0);
  } else {
    double integral
        = integrate(std::bind<double>(f, std::placeholders::_1, value_of(theta),
                                      x_r, x_i, std::ref(msgs)),
                    a, b, tolerance);

    size_t N = theta.size();
    std::vector<double> dintegral_dtheta(N);
    std::vector<double> theta_vals = value_of(theta);

    for (size_t n = 0; n < N; ++n) {
      dintegral_dtheta[n] = integrate(
          std::bind<double>(gradient_of_f<F>, f, std::placeholders::_1,
                            theta_vals, x_r, x_i, n, std::ref(msgs)),
          a, b, tolerance);
    }

    return precomputed_gradients(integral, theta, dintegral_dtheta);
  }
}

}  // namespace math
}  // namespace stan

#endif
