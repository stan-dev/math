#ifndef STAN_MATH_REV_FUNCTOR_integrate_1d_HPP
#define STAN_MATH_REV_FUNCTOR_integrate_1d_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/integrate_1d_adjoint.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/integrate_1d.hpp>
#include <cmath>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the integral of f from a to b to the given relative tolerance
 *
 * @tparam F Type of f
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam Args types of parameter pack arguments
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_var<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_impl(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    std::ostream *msgs, const Args &... args) {
  static constexpr const char *function = "integrate_1d";
  check_less_or_equal(function, "lower limit", a, b);
  return internal::integrate_1d_adjoint(
      function, f, a, b,
      [&](auto &&integrand) {
        return integrate(std::forward<decltype(integrand)>(integrand),
                         value_of(a), value_of(b), relative_tolerance);
      },
      msgs, args...);
}

/**
 * Compute the integral of the single variable function f from a to b to within
 * a specified relative tolerance. a and b can be finite or infinite.
 *
 * f should be compatible with reverse mode autodiff and have the signature:
 *   var f(double x, double xc, const std::vector<var>& theta,
 *     const std::vector<double>& x_r, const std::vector<int> &x_i,
 * std::ostream* msgs)
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream.
 *
 * Integrals that cross zero are broken into two, and the separate integrals are
 * each integrated to the given relative tolerance.
 *
 * For integrals with finite limits, the xc argument is the distance to the
 * nearest boundary. So for a > 0, b > 0, it will be a - x for x closer to a,
 * and b - x for x closer to b. xc is computed in a way that avoids the
 * precision loss of computing a - x or b - x manually. For integrals that cross
 * zero, xc can take values a - x, -x, or b - x depending on which integration
 * limit it is nearest.
 *
 * If either limit is infinite, xc is set to NaN
 *
 * The integration algorithm terminates when
 *   \f[
 *     \frac{{|I_{n + 1} - I_n|}}{{|I|_{n + 1}}} < \text{relative tolerance}
 *   \f]
 * where \f$I_{n}\f$ is the nth estimate of the integral and \f$|I|_{n}\f$ is
 * the nth estimate of the norm of the integral.
 *
 * Integrals that cross zero are
 * split into two. In this case, each integral is separately integrated to the
 * given relative_tolerance.
 *
 * Gradients of f that evaluate to NaN when the function evaluates to zero are
 * set to zero themselves. This is due to the autodiff easily overflowing to NaN
 * when evaluating gradients near the maximum and minimum floating point values
 * (where the function should be zero anyway for the integral to exist)
 *
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam T_theta type of parameters
 * @tparam T Type of f
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param theta additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename T_theta,
          typename = require_any_var_t<T_a, T_b, T_theta>>
inline return_type_t<T_a, T_b, T_theta> integrate_1d(
    const F &f, const T_a &a, const T_b &b, const std::vector<T_theta> &theta,
    const std::vector<double> &x_r, const std::vector<int> &x_i,
    std::ostream *msgs, const double relative_tolerance = std::sqrt(EPSILON)) {
  return integrate_1d_impl(integrate_1d_adapter<F>(f), a, b, relative_tolerance,
                           msgs, theta, x_r, x_i);
}

}  // namespace math
}  // namespace stan

#endif
