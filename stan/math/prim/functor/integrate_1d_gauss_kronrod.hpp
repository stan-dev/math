#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_GAUSS_KRONROD_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_1D_GAUSS_KRONROD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/integrate_1d_adapter.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <queue>
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
 * Integrate a single variable function f from a to b using Boost's
 * Gauss-Kronrod (G10,K21) rule and global adaptive error control. Subdivide
 * the interval with the largest estimated error until
 *   error <= max(relative_tolerance * L1, absolute_tolerance)
 * where error and L1 are the sums of the current panels' error and L1
 * estimates. Absolute tolerance participates in refinement, including for
 * derivative integrands. Failure to reach this target throws std::domain_error.
 *
 * The global subdivision strategy follows Algorithm 2 of P. Gonnet,
 * "A Review of Error Estimation in Adaptive Quadrature" (2010),
 * https://arxiv.org/abs/1003.4629. Mixed tolerances and the stagnation check
 * are adapted from QUADPACK's DQAG: https://www.netlib.org/quadpack/dqage.f.
 * This is not a full DQAG implementation: it retains Boost's embedded-rule
 * error estimate and Stan's L1 relative scale instead of DQAG's |Q| scale.
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
 * Infinite limits use the same rational changes of variable as Boost's
 * gauss_kronrod. The Jacobian scales the value, error, and L1 estimates.
 *
 * @tparam F Type of f
 * @param f the function to be integrated
 * @param a lower limit of integration (may be -infinity)
 * @param b upper limit of integration (may be +infinity)
 * @param relative_tolerance target error relative to the L1 estimate
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_depth maximum bisection depth of each panel
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate_gk(const F& f, double a, double b,
                           double relative_tolerance, double absolute_tolerance,
                           int max_depth) {
  static constexpr const char* function = "integrate_1d_gauss_kronrod";
  const bool infinite_a = std::isinf(a);
  const bool infinite_b = std::isinf(b);
  auto transformed_f = [&](double x) -> double {
    if (infinite_a && infinite_b) {
      const double inv = 1.0 / (1.0 - x * x);
      return f(x * inv, NOT_A_NUMBER) * inv * inv * (1.0 + x * x);
    } else if (infinite_a || infinite_b) {
      const double inv = 1.0 / (1.0 + x);
      const double offset = 2.0 * inv - 1.0;
      return f(infinite_a ? b - offset : a + offset, NOT_A_NUMBER) * inv * inv
             * 2.0;
    }
    return f(x, NOT_A_NUMBER);
  };
  struct panel {
    double a, b, value, error, l1;
    int depth;
  };
  auto evaluate = [&](double lower, double upper, int depth) {
    const double midpoint = lower / 2.0 + upper / 2.0;
    const double half_width = upper / 2.0 - lower / 2.0;
    auto mapped_f = [&](double x) {
      const double y = transformed_f(midpoint + half_width * x) * half_width;
      check_finite(function, "mapped integrand", y);
      return y;
    };
    panel result{lower, upper, 0.0, 0.0, 0.0, depth};
    // Integrate on [-1, 1] with depth zero. Including the Jacobian in
    // mapped_f also scales Boost 1.87's otherwise unscaled error estimate.
    using quadrature = boost::math::quadrature::gauss_kronrod<
        double, INTEGRATE_1D_GAUSS_KRONROD_ORDER>;
    result.value = quadrature::integrate(
        mapped_f, -1.0, 1.0, 0, relative_tolerance, &result.error, &result.l1);
    check_finite(function, "panel integral", result.value);
    check_finite(function, "panel error estimate", result.error);
    check_finite(function, "panel L1 estimate", result.l1);
    // QUADPACK's round-off floor accounts for summation even when the
    // embedded rules agree or the signed integral is close to zero.
    // See https://www.netlib.org/quadpack/dqk21.f.
    result.error = std::max(result.error, 50.0 * EPSILON * result.l1);
    return result;
  };
  const bool infinite = infinite_a || infinite_b;
  const panel first = evaluate(infinite ? -1.0 : a, infinite ? 1.0 : b, 0);
  double value = first.value;
  double error = first.error;
  double l1 = first.l1;
  auto converged = [&]() {
    return error <= std::max(relative_tolerance * l1, absolute_tolerance);
  };
  if (converged()) {
    return value;
  }

  // Keep every leaf so that incremental sums can be checked before returning.
  // Only leaves that can still be subdivided belong in the priority queue.
  std::vector<panel> panels{first};
  auto smaller_error = [&](std::size_t i, std::size_t j) {
    return panels[i].error < panels[j].error;
  };
  std::priority_queue<std::size_t, std::vector<std::size_t>,
                      decltype(smaller_error)>
      pending(smaller_error);
  if (max_depth > 0) {
    pending.push(0);
  }
  auto resum = [&]() {
    value = error = l1 = 0.0;
    for (const auto& leaf : panels) {
      value += leaf.value;
      error += leaf.error;
      l1 += leaf.l1;
    }
    check_finite(function, "integral", value);
    check_finite(function, "error estimate", error);
    check_finite(function, "L1 estimate", l1);
  };
  int roundoff_count = 0;
  while (true) {
    if (converged()) {
      resum();
      if (converged()) {
        return value;
      }
    }
    if (pending.empty() || roundoff_count >= 6
        || std::max(relative_tolerance * l1, absolute_tolerance)
               < 50.0 * EPSILON * l1) {
      resum();
      if (converged()) {
        return value;
      }
      throw_domain_error(
          function, "error estimate of integral", error, "",
          " exceeds max(relative_tolerance * L1, absolute_tolerance); "
          "bisection or floating-point precision limit reached");
    }
    const std::size_t index = pending.top();
    pending.pop();
    const panel parent = panels[index];
    const double midpoint = parent.a / 2.0 + parent.b / 2.0;
    if (midpoint == parent.a || midpoint == parent.b) {
      continue;
    }
    const panel left = evaluate(parent.a, midpoint, parent.depth + 1);
    const panel right = evaluate(midpoint, parent.b, parent.depth + 1);
    const double child_value = left.value + right.value;
    const double child_error = left.error + right.error;
    // DQAG's first round-off counter: six subdivisions with almost no
    // change in the integral and no reduction in estimated error.
    // Boost does not use DQAG's resasc-based error rescaling, so the
    // associated resasc saturation guard does not apply here.
    if (std::abs(parent.value - child_value) <= 1e-5 * std::abs(child_value)
        && child_error >= 0.99 * parent.error) {
      ++roundoff_count;
    }
    value += child_value - parent.value;
    error += child_error - parent.error;
    l1 += left.l1 + right.l1 - parent.l1;
    panels[index] = left;
    panels.push_back(right);
    if (left.depth < max_depth) {
      pending.push(index);
      pending.push(panels.size() - 1);
    }
  }
}

/**
 * Compute the integral of the single variable function f from a to b to within
 * a specified relative tolerance using adaptive Gauss-Kronrod (G10,K21)
 * quadrature. a and b can be finite or infinite.
 *
 * @tparam F type of function to integrate
 * @tparam Args types of additional arguments forwarded to f (all arithmetic)
 *
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance target error relative to the L1 estimate
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_depth maximum bisection depth of each panel
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments passed to f
 * @return numeric integral of function f
 */
template <typename F, typename... Args,
          require_all_st_arithmetic<Args...>* = nullptr>
inline double integrate_1d_gauss_kronrod_tol(const F& f, double a, double b,
                                             double relative_tolerance,
                                             double absolute_tolerance,
                                             int max_depth, std::ostream* msgs,
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
 * adaptive Gauss-Kronrod (G10,K21) quadrature. a and b can be finite or
 * infinite.
 *
 * The signature for f should be:
 *   double f(double x, double xc, std::ostream* msgs, Args... args...)
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream. xc is unused (always NaN) here; see
 * integrate_gk above for details.
 *
 * The integration algorithm terminates when the sum of the panels'
 * estimated quadrature errors satisfies
 *   \f[
 *     \text{error} \leq \max(\text{relative\_tolerance} \cdot |I|,
 *                            \text{absolute\_tolerance})
 *   \f]
 * where \f$|I|\f$ is the Boost estimate of the L1 norm of the integral.
 *
 *
 * @tparam F type of function to integrate
 * @tparam Args types of additional arguments forwarded to f (all arithmetic)
 *
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments passed to f
 * @return numeric integral of function f
 */

template <typename F, typename... Args,
          require_all_st_arithmetic<Args...>* = nullptr>
inline double integrate_1d_gauss_kronrod(const F& f, double a, double b,
                                         std::ostream* msgs,
                                         const Args&... args) {
  return integrate_1d_gauss_kronrod_tol(f, a, b, std::sqrt(EPSILON), 0.0,
                                        INTEGRATE_1D_GAUSS_KRONROD_MAX_DEPTH,
                                        msgs, args...);
}

}  // namespace math
}  // namespace stan

#endif
