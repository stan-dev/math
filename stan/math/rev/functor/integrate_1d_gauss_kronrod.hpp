#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_GAUSS_KRONROD_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_GAUSS_KRONROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/integrate_1d_adjoint.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/integrate_1d_gauss_kronrod.hpp>
#include <cmath>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the integral of f from a to b using adaptive Gauss-Kronrod (G21,K21)
 * quadrature.
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
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_depth maximum recursive bisection depth passed to Boost
 *   quadrature
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_var<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_gauss_kronrod_tol(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    double absolute_tolerance, int max_depth, std::ostream *msgs,
    const Args &... args) {
  static constexpr const char *function = "integrate_1d_gauss_kronrod";
  check_less_or_equal(function, "lower limit", a, b);
  check_nonnegative(function, "max_depth", max_depth);
  check_nonnegative(function, "absolute_tolerance", absolute_tolerance);
  return internal::integrate_1d_adjoint(
      function, f, a, b,
      [&](auto &&integrand) {
        return integrate_gk(std::forward<decltype(integrand)>(integrand),
                            value_of(a), value_of(b), relative_tolerance,
                            absolute_tolerance, max_depth);
      },
      msgs, args...);
}

/**
 * Compute the integral of the single variable function f from a to b using
 * adaptive Gauss-Kronrod (G21,K21) quadrature. a and b can be finite or
 * infinite.
 *
 * f should be compatible with reverse mode autodiff and have the signature:
 *   var f(double x, double xc, std::ostream* msgs, Args... args...);
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream. xc is unused (always NaN) here.
 *
 * The integration algorithm terminates when the Boost estimate of the
 * quadrature error satisfies
 *   error <= max(relative_tolerance * L1, absolute_tolerance)
 * where L1 is the Boost estimate of the L1 norm of the integral.
 *
 * Gradients of f that evaluate to NaN when the function evaluates to zero are
 * set to zero themselves. This is due to the autodiff easily overflowing to
 * NaN when evaluating gradients near the maximum and minimum floating point
 * values (where the function should be zero anyway for the integral to
 * exist).
 *
 * @tparam F Type of f
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam Args types of parameter pack arguments
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_var<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_gauss_kronrod(
    const F &f, const T_a &a, const T_b &b, std::ostream *msgs,
    const Args &... args) {
  return integrate_1d_gauss_kronrod_tol(f, a, b, std::sqrt(EPSILON), 0.0,
                                        INTEGRATE_1D_GAUSS_KRONROD_MAX_DEPTH,
                                        msgs, args...);
}

}  // namespace math
}  // namespace stan

#endif
