c#ifndef STAN_MATH_REV_ARR_FUNCTOR_integrate_1d_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_integrate_1d_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/arr/fun/to_var.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <cstdint>
#include <functional>
#include <cmath>
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
                            const std::vector<double>& x_r, const std::vector<int>& x_i,
                            size_t n, std::ostream& msgs) {
  double gradient = 0.0;
  start_nested();
  std::vector<var> theta_var(theta_vals.size());
  try {
    for(size_t i = 0; i < theta_vals.size(); i++)
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
 * Compute the integral of the single variable function f from a to b to within a
 * specified tolerance. a and b can be finite or infinite.
 *
 * a and b cannot be parameters.
 *
 * f should have the signature:
 *   double (double, std::vector<var>, std::vector<double>, std::vector<int>, std::ostream&)
 *
 * The first argument is the variable over which f is integrated.
 * The second argument are parameters, which should be a vector of stan::math::var
 * The third argument is extra real data, which should always be a vector of doubles
 * The fourth argument is extra integer data, which should always be a vector of ints
 * The fifth argument is a stream to which error messages are printed
 *
 * Boost decides the integration is converged when subsequent estimates of the
 * integral are less than tolerance * abs(integral). This means the tolerance is
 * relative to the actual scale of the integral.
 *
 * @tparam T Type of f
 * @param f a functor with signature
 *   double (double, std::vector<var>, std::vector<double>, std::vector<int>, std::ostream&)
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param param additional parameters to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param tolerance integrator tolerance passed to Boost quadrature
 * @return numeric integral of function f
 */
template <typename F>
inline var
integrate_1d(const F& f, const double a, const double b, const std::vector<var>& theta,
             const std::vector<double>& x_r, const std::vector<int>& x_i,
    std::ostream& msgs, const double tolerance = 1e-6) {
  static const char *function = "integrate_1d";
  check_less_or_equal(function, "lower limit", a, b);

  if(a == b) {
    return var(0.0);
  } else {
    double integral =
        integrate(std::bind<double>(f, std::placeholders::_1, value_of(theta), x_r, x_i, std::ref(msgs)),
                  a, b, tolerance);

    size_t N = theta.size();
    std::vector<double> dintegral_dtheta(N);
    std::vector<double> theta_vals = value_of(theta);

    for (size_t n = 0; n < N; ++n) {
      dintegral_dtheta[n] = integrate(std::bind<double>(gradient_of_f<F>, f, std::placeholders::_1, theta_vals, x_r, x_i, n, std::ref(msgs)), a, b, tolerance);
    }

    return precomputed_gradients(integral, theta, dintegral_dtheta);
  }
}

}  // namespace math
}  // namespace stan

#endif
