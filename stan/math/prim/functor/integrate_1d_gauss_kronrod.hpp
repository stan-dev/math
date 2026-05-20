#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_GAUSS_KRONROD_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_GAUSS_KRONROD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/integrate_1d_adapter.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Default Kronrod order used by integrate_1d_gauss_kronrod. Boost provides
 * compile-time tables for N in {15, 21, 31, 41, 51, 61}; 21 is the common
 * QUADPACK choice and a reasonable speed/accuracy trade-off for smooth
 * integrands.
 */
constexpr unsigned int INTEGRATE_1D_GAUSS_KRONROD_ORDER = 21;

/**
 * Default recursive bisection depth used by integrate_1d_gauss_kronrod
 * when the user does not pass one explicitly. Matches Boost's default.
 */
constexpr int INTEGRATE_1D_GAUSS_KRONROD_MAX_DEPTH = 15;

/**
 * Integrate a single variable function f from a to b using Boost's adaptive
 * Gauss-Kronrod (G21,K21) quadrature, with QUADPACK-style mixed convergence
 * criterion. The integration succeeds (returns the Boost estimate Q)
 * whenever
 *   error <= max(relative_tolerance * L1, absolute_tolerance)
 * where error and L1 are Boost's quadrature-error and L1-norm estimates.
 * A larger error throws std::domain_error.
 *
 * Setting absolute_tolerance to a small positive value lets callers escape
 * the pathological regime where the relative-tolerance test on its own is
 * checking accumulated floating-point round-off against itself (this
 * happens routinely in nested integrate_1d_gauss_kronrod calls when the
 * outer integration probes the deep tail of the integrand and every
 * inner evaluation sees an essentially-zero integrand). Setting it to
 * zero (the default) reproduces the strict pure-relative-tolerance
 * behaviour of integrate_1d.
 *
 * The signature for f should be:
 *   double f(double x, double xc)
 *
 * Unlike integrate_1d (which uses tanh_sinh/exp_sinh/sinh_sinh and computes a
 * meaningful distance-to-boundary xc), Gauss-Kronrod does not produce xc, so
 * this routine always passes xc == NaN to the user functor. User functors
 * written for integrate_1d that rely on xc must be rewritten without it before
 * being used here.
 *
 * Boost's gauss_kronrod handles infinite limits internally via the usual
 * change of variable; no special handling for integrals crossing zero is
 * required.
 *
 * @tparam F Type of f
 * @param f the function to be integrated
 * @param a lower limit of integration (may be -infinity)
 * @param b upper limit of integration (may be +infinity)
 * @param relative_tolerance target relative tolerance passed to Boost
 * quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_depth maximum recursive bisection depth passed to Boost
 * quadrature
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate_gk(const F& f, double a, double b,
                           double relative_tolerance,
                           double absolute_tolerance, int max_depth) {
  static constexpr const char* function = "integrate_1d_gauss_kronrod";
  double error = 0.0;
  double L1 = 0.0;

  // Gauss-Kronrod does not pass a distance-to-boundary; the user functor
  // still takes (x, xc) for signature compatibility with integrate_1d, but
  // xc is unused here.
  auto f_wrap = [&f](double x) { return f(x, NOT_A_NUMBER); };

  using boost::math::quadrature::gauss_kronrod;
  const unsigned int depth = max_depth < 0
                                 ? 0u
                                 : static_cast<unsigned int>(max_depth);
  double Q = gauss_kronrod<double, INTEGRATE_1D_GAUSS_KRONROD_ORDER>::integrate(
      f_wrap, a, b, depth, relative_tolerance, &error, &L1);

  // QUADPACK-style mixed convergence: throw only if the Boost error
  // exceeds both the relative-tolerance target (rel_tol * L1) and the
  // user-supplied absolute floor. With absolute_tolerance = 0 (default)
  // this reduces to the strict pure-relative test.
  const double convergence_threshold
      = std::max(relative_tolerance * L1, absolute_tolerance);
  if (error > convergence_threshold) {
    [error]() STAN_COLD_PATH {
      throw_domain_error(
          function, "error estimate of integral", error, "",
          " exceeds max(relative_tolerance * L1, absolute_tolerance)");
    }();
  }
  return Q;
}

/**
 * Compute the integral of the single variable function f from a to b to within
 * a specified relative tolerance using adaptive Gauss-Kronrod (G21,K21)
 * quadrature. a and b can be finite or infinite.
 *
 * @tparam F type of function to integrate
 * @tparam Args types of additional arguments forwarded to f (all arithmetic)
 *
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_depth maximum recursive bisection depth passed to Boost quadrature
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments passed to f
 * @return numeric integral of function f
 */
template <typename F, typename... Args,
          require_all_st_arithmetic<Args...>* = nullptr>
inline double integrate_1d_gauss_kronrod_impl(
    const F& f, double a, double b, double relative_tolerance,
    double absolute_tolerance, int max_depth, std::ostream* msgs,
    const Args&... args) {
  static constexpr const char* function = "integrate_1d_gauss_kronrod";
  check_less_or_equal(function, "lower limit", a, b);
  check_nonnegative(function, "max_depth", max_depth);
  check_nonnegative(function, "absolute_tolerance", absolute_tolerance);
  if (unlikely(a == b)) {
    if (std::isinf(a)) {
      throw_domain_error(function, "Integration endpoints are both", a, "", "");
    }
    return 0.0;
  } else {
    return integrate_gk(
        [&](auto&& x, auto&& xc) { return f(x, xc, msgs, args...); }, a, b,
        relative_tolerance, absolute_tolerance, max_depth);
  }
}

/**
 * Compute the integral of the single variable function f from a to b using
 * adaptive Gauss-Kronrod (G21,K21) quadrature. a and b can be finite or
 * infinite.
 *
 * The signature for f should be:
 *   double f(double x, double xc, const std::vector<double>& theta,
 *     const std::vector<double>& x_r, const std::vector<int>& x_i,
 *     std::ostream* msgs)
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream. xc is unused (always NaN) here; see
 * integrate_gk above for details.
 *
 * The integration algorithm terminates when the Boost estimate of the
 * quadrature error satisfies
 *   \f[
 *     \text{error} \leq \max(\text{relative\_tolerance} \cdot |I|,
 *                            \text{absolute\_tolerance})
 *   \f]
 * where \f$|I|\f$ is the Boost estimate of the L1 norm of the integral.
 *
 * @tparam F type of function to integrate
 *
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param theta additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_depth maximum recursive bisection depth passed to Boost
 *   quadrature
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate_1d_gauss_kronrod(
    const F& f, double a, double b, const std::vector<double>& theta,
    const std::vector<double>& x_r, const std::vector<int>& x_i,
    std::ostream* msgs,
    const double relative_tolerance = std::sqrt(EPSILON),
    const double absolute_tolerance = 0.0,
    const int max_depth = INTEGRATE_1D_GAUSS_KRONROD_MAX_DEPTH) {
  return integrate_1d_gauss_kronrod_impl(integrate_1d_adapter<F>(f), a, b,
                                         relative_tolerance,
                                         absolute_tolerance, max_depth, msgs,
                                         theta, x_r, x_i);
}

}  // namespace math
}  // namespace stan

#endif
