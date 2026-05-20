#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/integrate_1d_adapter.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Default maximum refinement count used by integrate_1d_double_exponential
 * when the user does not pass one explicitly. Matches Boost's tanh_sinh
 * default (Boost's exp_sinh/sinh_sinh default is 9; we use 15 here for
 * symmetry with integrate_1d_gauss_kronrod's max_depth = 15 and because
 * tanh_sinh is by far the most-dispatched branch).
 */
constexpr int INTEGRATE_1D_DOUBLE_EXPONENTIAL_MAX_REFINEMENTS = 15;

/**
 * Integrate a single variable function f from a to b using Boost's adaptive
 * double-exponential quadrature, with QUADPACK-style mixed convergence
 * criterion. The integration succeeds (returns the Boost estimate Q)
 * whenever
 *   error <= max(relative_tolerance * L1, absolute_tolerance)
 * on each piece, where error and L1 are Boost's quadrature-error and
 * L1-norm estimates for that piece. A larger error throws
 * std::domain_error.
 *
 * Setting absolute_tolerance to a small positive value escapes the
 * pathological regime where the relative-tolerance test on its own is
 * checking accumulated floating-point round-off against itself. This
 * happens in particular under the zero-crossing split (see below) when
 * one half of the split integrates to near-zero: the strict per-piece
 * relative test on that half can fire spuriously. Setting it to zero
 * (the default) reproduces the strict pure-relative-tolerance
 * behaviour of integrate_1d.
 *
 * The signature for f should be:
 *   double f(double x, double xc)
 *
 * Unlike integrate_1d_gauss_kronrod, double-exponential quadrature
 * computes a meaningful distance-to-boundary xc, and user functors may
 * (and should, for accuracy near singular endpoints) use it. For
 * a > 0, b > 0, xc is a - x for x closer to a, and b - x for x closer
 * to b. xc is computed in a way that avoids the precision loss of
 * computing a - x or b - x manually. For integrals that cross zero, xc
 * can take values a - x, -x, or b - x depending on which integration
 * limit is nearest. If either limit is infinite, xc is set to NaN.
 *
 * Depending on whether or not a is finite or negative infinity and b is
 * finite or positive infinity, a different version of the 1d quadrature
 * algorithm from the Boost quadrature library is chosen.
 *
 * Integrals that cross zero are broken into two, and the separate
 * integrals are each integrated to the given tolerances. This is the
 * pre-existing behaviour of integrate_1d, preserved here to maintain
 * call-site compatibility.
 *
 * @tparam F Type of f
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance target relative tolerance passed to Boost
 *   quadrature
 * @param absolute_tolerance absolute-error floor on the per-piece
 *   convergence test
 * @param max_refinements maximum refinement level passed to the
 *   constructor of the Boost quadrature class
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate_de(const F& f, double a, double b,
                           double relative_tolerance,
                           double absolute_tolerance, int max_refinements) {
  static constexpr const char* function = "integrate_1d_double_exponential";
  double error1 = 0.0;
  double error2 = 0.0;
  double L1 = 0.0;
  double L2 = 0.0;
  size_t levels = 0;

  const size_t mr = max_refinements < 0
                        ? 0u
                        : static_cast<size_t>(max_refinements);

  auto one_integral_convergence_check = [&](double e, double rel_tol,
                                            double abs_tol, double L) {
    const double threshold = std::max(rel_tol * L, abs_tol);
    if (e > threshold) {
      [e]() STAN_COLD_PATH {
        throw_domain_error(
            function, "error estimate of integral", e, "",
            " exceeds max(relative_tolerance * L1, absolute_tolerance)");
      }();
    }
  };
  auto two_integral_convergence_check = [&](double e1, double e2,
                                            double rel_tol, double abs_tol,
                                            double La, double Lb) {
    const double threshold_a = std::max(rel_tol * La, abs_tol);
    const double threshold_b = std::max(rel_tol * Lb, abs_tol);
    if (e1 > threshold_a) {
      [e1]() STAN_COLD_PATH {
        throw_domain_error(
            function, "error estimate of integral below zero", e1, "",
            " exceeds max(relative_tolerance * L1, absolute_tolerance) for "
            "the lower piece of the split");
      }();
    }
    if (e2 > threshold_b) {
      [e2]() STAN_COLD_PATH {
        throw_domain_error(
            function, "error estimate of integral above zero", e2, "",
            " exceeds max(relative_tolerance * L1, absolute_tolerance) for "
            "the upper piece of the split");
      }();
    }
  };

  // If a or b is infinite, set xc to NaN (see docs above)
  auto f_wrap = [&f](double x) { return f(x, NOT_A_NUMBER); };
  if (std::isinf(a) && std::isinf(b)) {
    boost::math::quadrature::sinh_sinh<double> integrator(mr);
    double Q = integrator.integrate(f_wrap, relative_tolerance, &error1, &L1,
                                    &levels);
    one_integral_convergence_check(error1, relative_tolerance,
                                   absolute_tolerance, L1);
    return Q;
  } else if (std::isinf(a)) {
    boost::math::quadrature::exp_sinh<double> integrator(mr);
    // If the integral crosses zero, break it into two (advice from the
    // Boost implementation).
    if (b <= 0.0) {
      double Q = integrator.integrate(f_wrap, a, b, relative_tolerance, &error1,
                                      &L1, &levels);
      one_integral_convergence_check(error1, relative_tolerance,
                                     absolute_tolerance, L1);
      return Q;
    } else {
      boost::math::quadrature::tanh_sinh<double> integrator_right(mr);
      double Q = integrator.integrate(f_wrap, a, 0.0, relative_tolerance,
                                      &error1, &L1, &levels)
                 + integrator_right.integrate(
                     f_wrap, 0.0, b, relative_tolerance, &error2, &L2, &levels);
      two_integral_convergence_check(error1, error2, relative_tolerance,
                                     absolute_tolerance, L1, L2);
      return Q;
    }
  } else if (std::isinf(b)) {
    boost::math::quadrature::exp_sinh<double> integrator(mr);
    if (a >= 0.0) {
      double Q = integrator.integrate(f_wrap, a, b, relative_tolerance, &error1,
                                      &L1, &levels);
      one_integral_convergence_check(error1, relative_tolerance,
                                     absolute_tolerance, L1);
      return Q;
    } else {
      boost::math::quadrature::tanh_sinh<double> integrator_left(mr);
      double Q = integrator_left.integrate(f_wrap, a, 0, relative_tolerance,
                                           &error1, &L1, &levels)
                 + integrator.integrate(f_wrap, relative_tolerance, &error2,
                                        &L2, &levels);
      two_integral_convergence_check(error1, error2, relative_tolerance,
                                     absolute_tolerance, L1, L2);
      return Q;
    }
  } else {
    auto f_wrap = [&f](double x, double xc) { return f(x, xc); };
    boost::math::quadrature::tanh_sinh<double> integrator(mr);
    if (a < 0.0 && b > 0.0) {
      double Q = integrator.integrate(f_wrap, a, 0.0, relative_tolerance,
                                      &error1, &L1, &levels)
                 + integrator.integrate(f_wrap, 0.0, b, relative_tolerance,
                                        &error2, &L2, &levels);
      two_integral_convergence_check(error1, error2, relative_tolerance,
                                     absolute_tolerance, L1, L2);
      return Q;
    } else {
      double Q = integrator.integrate(f_wrap, a, b, relative_tolerance, &error1,
                                      &L1, &levels);
      one_integral_convergence_check(error1, relative_tolerance,
                                     absolute_tolerance, L1);
      return Q;
    }
  }
}

/**
 * Compute the integral of the single variable function f from a to b using
 * adaptive double-exponential quadrature. a and b can be finite or infinite.
 *
 * @tparam F type of function to integrate
 * @tparam Args types of additional arguments forwarded to f (all arithmetic)
 *
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_refinements maximum refinement level passed to the
 *   Boost quadrature class constructor
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments passed to f
 * @return numeric integral of function f
 */
template <typename F, typename... Args,
          require_all_st_arithmetic<Args...>* = nullptr>
inline double integrate_1d_double_exponential_impl(
    const F& f, double a, double b, double relative_tolerance,
    double absolute_tolerance, int max_refinements, std::ostream* msgs,
    const Args&... args) {
  static constexpr const char* function = "integrate_1d_double_exponential";
  check_less_or_equal(function, "lower limit", a, b);
  check_nonnegative(function, "max_refinements", max_refinements);
  check_nonnegative(function, "absolute_tolerance", absolute_tolerance);
  if (unlikely(a == b)) {
    if (std::isinf(a)) {
      throw_domain_error(function, "Integration endpoints are both", a, "", "");
    }
    return 0.0;
  } else {
    return integrate_de(
        [&](auto&& x, auto&& xc) { return f(x, xc, msgs, args...); }, a, b,
        relative_tolerance, absolute_tolerance, max_refinements);
  }
}

/**
 * Compute the integral of the single variable function f from a to b using
 * adaptive double-exponential quadrature. a and b can be finite or infinite.
 *
 * The signature for f should be:
 *   double f(double x, double xc, const std::vector<double>& theta,
 *     const std::vector<double>& x_r, const std::vector<int>& x_i,
 *     std::ostream* msgs)
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream.
 *
 * The integration algorithm terminates per piece when
 *   error <= max(relative_tolerance * L1, absolute_tolerance)
 * where L1 is the Boost estimate of the L1 norm of the integral.
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
 * @param relative_tolerance tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_refinements maximum refinement level passed to the Boost
 *   quadrature class constructor
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate_1d_double_exponential(
    const F& f, double a, double b, const std::vector<double>& theta,
    const std::vector<double>& x_r, const std::vector<int>& x_i,
    std::ostream* msgs,
    const double relative_tolerance = std::sqrt(EPSILON),
    const double absolute_tolerance = 0.0,
    const int max_refinements
    = INTEGRATE_1D_DOUBLE_EXPONENTIAL_MAX_REFINEMENTS) {
  return integrate_1d_double_exponential_impl(
      integrate_1d_adapter<F>(f), a, b, relative_tolerance, absolute_tolerance,
      max_refinements, msgs, theta, x_r, x_i);
}

}  // namespace math
}  // namespace stan

#endif
