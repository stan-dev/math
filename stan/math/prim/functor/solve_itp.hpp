#ifndef STAN_MATH_PRIM_FUNCTOR_SOLVE_ITP_HPP
#define STAN_MATH_PRIM_FUNCTOR_SOLVE_ITP_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace stan {
namespace math {
namespace internal {

/**
 * Return ceil(log2(width / (2 * x_tolerance))), truncated below at zero,
 * without forming a potentially overflowing or underflowing ratio.
 */
inline int64_t itp_bisection_iterations(double width, double x_tolerance) {
  if (width <= 2.0 * x_tolerance) {
    return 0;
  }

  int width_exponent;
  int tolerance_exponent;
  const double width_mantissa = std::frexp(width, &width_exponent);
  const double tolerance_mantissa
      = std::frexp(x_tolerance, &tolerance_exponent);

  // ceil(log2(width / x_tolerance)) is the exponent difference plus one
  // exactly when the mantissa ratio is greater than one.  Subtracting one
  // then gives ceil(log2(width / (2 * x_tolerance))).
  const int64_t ceil_log2_ratio
      = static_cast<int64_t>(width_exponent)
        - static_cast<int64_t>(tolerance_exponent)
        + (width_mantissa > tolerance_mantissa ? 1 : 0);
  return std::max<int64_t>(0, ceil_log2_ratio - 1);
}

/**
 * Scale x_tolerance by an integral power of two.  The result is allowed to
 * saturate at positive infinity because it is subsequently used only as an
 * upper projection radius.
 */
inline double itp_scaled_tolerance(double x_tolerance, int64_t exponent) {
  if (exponent > std::numeric_limits<int>::max()) {
    return std::numeric_limits<double>::infinity();
  }
  if (exponent < std::numeric_limits<int>::min()) {
    return 0.0;
  }
  return std::ldexp(x_tolerance, static_cast<int>(exponent));
}

/**
 * Return the default truncation scale without overflowing on an interval whose
 * positive width is subnormal.  The fallback value is irrelevant when the
 * default tolerance already contains the entire interval, but keeping it finite
 * allows endpoint and bracket checks to retain their normal semantics.
 */
inline double itp_default_kappa1(double width) {
  const double overflow_threshold = 0.2 / std::numeric_limits<double>::max();
  return width > overflow_threshold ? 0.2 / width : 1.0;
}

/**
 * Return the regula-falsi point for an opposite-signed bracket.  This form is
 * algebraically equivalent to
 *
 *   (f_b * a - f_a * b) / (f_b - f_a),
 *
 * but avoids products of a coordinate and a potentially extreme function
 * value, and avoids overflow in |f_a| + |f_b|.
 */
inline double itp_false_position(double a, double b, double f_a, double f_b) {
  const double abs_f_a = std::abs(f_a);
  const double abs_f_b = std::abs(f_b);
  double lambda;

  if (abs_f_a <= abs_f_b) {
    const double ratio = abs_f_a / abs_f_b;
    lambda = ratio / (1.0 + ratio);
  } else {
    const double ratio = abs_f_b / abs_f_a;
    lambda = 1.0 / (1.0 + ratio);
  }
  return a + lambda * (b - a);
}

inline bool itp_opposite_signs(double x, double y) {
  return (x < 0.0 && y > 0.0) || (x > 0.0 && y < 0.0);
}

inline void itp_check_controls(const char* function, double lower, double upper,
                               double x_tolerance, double kappa1, double kappa2,
                               int64_t n0, int64_t max_num_steps) {
  static constexpr double one_plus_phi
      = 0.5 * (3.0 + 2.23606797749978969640917366873127623544);

  check_finite(function, "lower bound", lower);
  check_finite(function, "upper bound", upper);
  check_less(function, "lower bound", lower, upper);
  check_finite(function, "x tolerance", x_tolerance);
  check_positive(function, "x tolerance", x_tolerance);
  check_finite(function, "kappa1", kappa1);
  check_positive(function, "kappa1", kappa1);
  check_finite(function, "kappa2", kappa2);
  check_greater_or_equal(function, "kappa2", kappa2, 1.0);
  check_less(function, "kappa2", kappa2, one_plus_phi);
  check_nonnegative(function, "n0", n0);
  check_positive(function, "max_num_steps", max_num_steps);

  const double width = upper - lower;
  check_finite(function, "upper bound minus lower bound", width);
}

/**
 * Primitive ITP implementation after validation of control arguments.
 */
template <typename F, typename... Args>
inline double solve_itp_primitive(const F& f, double lower, double upper,
                                  double x_tolerance, double kappa1,
                                  double kappa2, int64_t n0,
                                  int64_t max_num_steps, std::ostream* msgs,
                                  const Args&... args) {
  static constexpr const char* function = "solve_itp";
  double a = lower;
  double b = upper;
  double f_a = value_of(f(a, msgs, args...));
  double f_b = value_of(f(b, msgs, args...));

  check_finite(function, "function value at lower bound", f_a);
  check_finite(function, "function value at upper bound", f_b);

  if (f_a == 0.0) {
    return a;
  }
  if (f_b == 0.0) {
    return b;
  }
  if (!itp_opposite_signs(f_a, f_b)) {
    std::ostringstream message;
    message << function
            << ": function values at the lower and upper bounds are " << f_a
            << " and " << f_b << ", but must have opposite signs";
    throw std::domain_error(message.str());
  }

  const int64_t n_half = itp_bisection_iterations(b - a, x_tolerance);
  if (n0 > std::numeric_limits<int64_t>::max() - n_half) {
    throw_domain_error(function, "n0", n0, "",
                       " is too large to form the ITP iteration bound");
  }
  const int64_t n_max = n_half + n0;
  if (n_max > max_num_steps) {
    throw_domain_error(function, "theoretical ITP iteration bound", n_max, "",
                       " exceeds max_num_steps");
  }

  int64_t iteration = 0;
  const double two_tolerance = 2.0 * x_tolerance;
  while (b - a > two_tolerance) {
    // The real-arithmetic proof gives n_max iterations.  One additional
    // midpoint step is permitted solely to repair an outward-rounded final
    // width such as nextafter(2 * x_tolerance, +infinity).
    if (iteration > n_max) {
      throw_domain_error(
          function, "number of ITP iterations", iteration, "",
          " exceeded the exact-arithmetic bound plus one roundoff repair");
    }
    if (iteration >= max_num_steps) {
      throw_domain_error(function, "number of ITP iterations", iteration, "",
                         " reached max_num_steps before convergence");
    }

    const double width = b - a;
    const double midpoint = 0.5 * a + 0.5 * b;
    if (!(midpoint > a && midpoint < b)) {
      throw_domain_error(
          function, "x tolerance", x_tolerance, "",
          " is smaller than the representable spacing in the final bracket");
    }

    double x_itp = midpoint;
    if (iteration < n_max) {
      // Interpolation.
      const double x_f = itp_false_position(a, b, f_a, f_b);

      // Truncation toward the midpoint.
      const double gap = std::abs(midpoint - x_f);
      const double sigma = midpoint > x_f ? 1.0 : (midpoint < x_f ? -1.0 : 0.0);
      double x_t = midpoint;
      if (gap > 0.0) {
        const double log_delta = std::log(kappa1) + kappa2 * std::log(width);
        if (log_delta < std::log(gap)) {
          x_t = x_f + sigma * std::exp(log_delta);
        }
      }

      // Projection onto the minmax-safe interval centered at the midpoint.
      const double scaled_tolerance
          = itp_scaled_tolerance(x_tolerance, n_max - iteration);
      double projection_radius = scaled_tolerance - 0.5 * width;
      projection_radius = std::max(0.0, projection_radius);
      const double projected_distance
          = std::min(projection_radius, std::abs(x_t - midpoint));
      x_itp = midpoint - sigma * projected_distance;

      // Midpoint fallback protects the bracket from exceptional roundoff.
      if (!std::isfinite(x_itp) || !(x_itp > a && x_itp < b)) {
        x_itp = midpoint;
      }
    }

    const double f_itp = value_of(f(x_itp, msgs, args...));
    check_finite(function, "function value at candidate point", f_itp);
    if (f_itp == 0.0) {
      return x_itp;
    }

    if (itp_opposite_signs(f_a, f_itp)) {
      b = x_itp;
      f_b = f_itp;
    } else {
      a = x_itp;
      f_a = f_itp;
    }
    ++iteration;
  }

  return 0.5 * a + 0.5 * b;
}

}  // namespace internal

/**
 * Find a root of a continuous scalar function in a finite sign-changing
 * bracket using the interpolate-truncate-project (ITP) method.
 *
 * The callback must be callable as
 *
 *   f(x, msgs, args...)
 *
 * and return a scalar.  The interval endpoints and tuning controls are data;
 * differentiation with respect to args is supplied by the reverse-mode
 * overload.
 *
 * @tparam F callback type
 * @tparam Args types of additional callback arguments
 * @param f scalar callback
 * @param lower lower endpoint of the initial bracket
 * @param upper upper endpoint of the initial bracket
 * @param x_tolerance requested absolute error bound on the returned midpoint
 * @param kappa1 positive truncation scale
 * @param kappa2 truncation exponent in [1, 1 + phi)
 * @param n0 nonnegative slack in the bisection minmax iteration bound
 * @param max_num_steps maximum permitted number of function-evaluation
 *        iterations; one step beyond the exact-arithmetic ITP bound may be
 *        needed to repair an outward-rounded final width
 * @param msgs stream for callback messages
 * @param args additional callback arguments
 * @return midpoint of a final sign-changing bracket whose computed width is
 *         at most 2 * x_tolerance, or an exactly evaluated root
 */
template <typename F, typename... Args,
          require_all_st_arithmetic<Args...>* = nullptr>
inline double solve_itp_tol(const F& f, double lower, double upper,
                            double x_tolerance, double kappa1, double kappa2,
                            int64_t n0, int64_t max_num_steps,
                            std::ostream* msgs, const Args&... args) {
  static constexpr const char* function = "solve_itp";
  internal::itp_check_controls(function, lower, upper, x_tolerance, kappa1,
                               kappa2, n0, max_num_steps);
  return internal::solve_itp_primitive(f, lower, upper, x_tolerance, kappa1,
                                       kappa2, n0, max_num_steps, msgs,
                                       args...);
}

/**
 * Find a root using the default ITP controls
 *
 *   x_tolerance = 1e-10,
 *   kappa1       = 0.2 / (upper - lower),
 *   kappa2       = 2,
 *   n0           = 1,
 *   max steps    = 200.
 */
template <typename F, typename... Args,
          require_all_st_arithmetic<Args...>* = nullptr>
inline double solve_itp(const F& f, double lower, double upper,
                        std::ostream* msgs, const Args&... args) {
  static constexpr const char* function = "solve_itp";
  check_finite(function, "lower bound", lower);
  check_finite(function, "upper bound", upper);
  check_less(function, "lower bound", lower, upper);
  const double width = upper - lower;
  check_finite(function, "upper bound minus lower bound", width);

  return solve_itp_tol(f, lower, upper, 1e-10,
                       internal::itp_default_kappa1(width), 2.0, 1, 200, msgs,
                       args...);
}

}  // namespace math
}  // namespace stan

#endif
