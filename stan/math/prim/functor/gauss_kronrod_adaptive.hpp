#ifndef STAN_MATH_PRIM_FUNCTOR_GAUSS_KRONROD_ADAPTIVE_HPP
#define STAN_MATH_PRIM_FUNCTOR_GAUSS_KRONROD_ADAPTIVE_HPP

// The quadrature driver below is a transcription of Boost.Math's
// boost::math::quadrature::gauss_kronrod<Real, N>::integrate and its private
// helpers integrate_non_adaptive_m1_1 and recursive_adaptive_integrate
// (Copyright John Maddock 2017, Copyright Nick Thompson 2017), specialised to
// Real = double and N = 21 and used under the Boost Software License,
// Version 1.0.  See licenses/boost-license.txt.
//
// It exists because Boost's public integrate() hard-codes the recursion's
// absolute-tolerance budget to zero, and that budget is not otherwise
// reachable.  Everything else is reproduced verbatim so that
// absolute_tolerance == 0 is bit-for-bit identical to calling Boost; see
// gauss_kronrod_21_integrate below.

#include <stan/math/prim/err/throw_domain_error.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/precision.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace stan {
namespace math {
namespace internal {

/**
 * Apply the (G10, K21) pair to f over [-1, 1].
 *
 * Specialisation of Boost's integrate_non_adaptive_m1_1 to N = 21.  The Gauss
 * order (N - 1) / 2 = 10 is even, so the Gauss rule has no node at the origin:
 * the midpoint contributes to the Kronrod sum only, and the two loops start at
 * 1 and 2 respectively.
 *
 * The error estimate is |K21 - G10|, floored at the round-off level of the
 * Kronrod sum.  Following Boost, it is NOT rescaled by the panel half-width;
 * callers therefore receive an error in the units of the [-1, 1] reference
 * panel while the estimate and L1 norm are in the units of the original
 * interval.  This is preserved deliberately: integrate_1d_gauss_kronrod
 * already compares this quantity against absolute_tolerance, so changing the
 * convention would silently move a released throw threshold.
 *
 * @tparam F type of f
 * @param f function to integrate over [-1, 1]
 * @param[out] error estimate of the quadrature error (must not be null)
 * @param[out] l1 estimate of the L1 norm of f (must not be null)
 * @return the K21 estimate of the integral of f over [-1, 1]
 */
template <typename F>
inline double gauss_kronrod_21_m1_1(const F& f, double* error, double* l1) {
  using boost::math::quadrature::gauss;
  using boost::math::quadrature::gauss_kronrod;
  const auto& abscissa = gauss_kronrod<double, 21>::abscissa();
  const auto& kronrod_weights = gauss_kronrod<double, 21>::weights();
  const auto& gauss_weights = gauss<double, 10>::weights();

  double fp = f(0.0);
  double fm = 0.0;
  double kronrod_result = fp * static_cast<double>(kronrod_weights[0]);
  double gauss_result = 0.0;
  double L1 = std::abs(kronrod_result);

  for (unsigned int i = 1; i < abscissa.size(); i += 2) {
    fp = f(static_cast<double>(abscissa[i]));
    fm = f(static_cast<double>(-abscissa[i]));
    kronrod_result += (fp + fm) * static_cast<double>(kronrod_weights[i]);
    L1 += (std::abs(fp) + std::abs(fm))
          * static_cast<double>(kronrod_weights[i]);
    gauss_result += (fp + fm) * static_cast<double>(gauss_weights[i / 2]);
  }
  for (unsigned int i = 2; i < abscissa.size(); i += 2) {
    fp = f(static_cast<double>(abscissa[i]));
    fm = f(static_cast<double>(-abscissa[i]));
    kronrod_result += (fp + fm) * static_cast<double>(kronrod_weights[i]);
    L1 += (std::abs(fp) + std::abs(fm))
          * static_cast<double>(kronrod_weights[i]);
  }

  *l1 = L1;
  *error = std::max(
      std::abs(kronrod_result - gauss_result),
      std::abs(kronrod_result * boost::math::tools::epsilon<double>() * 2.0));
  return kronrod_result;
}

/**
 * Adaptively integrate f over the finite, ordered interval [a, b] by local
 * bisection.
 *
 * A panel is bisected only when its own error estimate exceeds BOTH a local
 * relative target (|estimate| * relative_tolerance) and the absolute budget
 * absolute_tolerance, and the budget is halved on the way down so that the
 * per-panel floors sum to the caller's floor.  Refinement is therefore driven
 * by each panel in isolation: work is bounded by the panels that actually need
 * it, and a panel that already meets either target is never touched again.
 *
 * At the root panel the budget is max(absolute_tolerance, |estimate| *
 * relative_tolerance).  Boost instead treats a zero budget as a sentinel and
 * overwrites it with the root's local relative target; taking the max
 * subsumes that (absolute_tolerance == 0 gives exactly Boost's value, so the
 * whole routine stays bit-for-bit Boost) while making the floor continuous in
 * absolute_tolerance: a negligible floor is a no-op rather than a request for
 * unbounded refinement.  The resulting rule, max(relative target, absolute
 * floor), is the same mixed QUADPACK criterion that integrate_1d_gauss_kronrod
 * applies to decide convergence, so refinement and acceptance agree.
 *
 * @tparam F type of f
 * @param f function to integrate
 * @param relative_tolerance local relative target per panel
 * @param a lower limit (finite, a < b)
 * @param b upper limit (finite)
 * @param max_levels remaining bisection levels
 * @param absolute_tolerance absolute error budget for this panel
 * @param is_root whether this is the top-level panel
 * @param[out] error estimate of the quadrature error (must not be null)
 * @param[out] l1 estimate of the L1 norm (must not be null)
 * @return estimate of the integral of f over [a, b]
 */
template <typename F>
inline double gauss_kronrod_21_recursive(const F& f, double relative_tolerance,
                                         double a, double b,
                                         unsigned int max_levels,
                                         double absolute_tolerance,
                                         bool is_root, double* error,
                                         double* l1) {
  double error_local = 0.0;
  const double mean = (b + a) / 2;
  const double scale = (b - a) / 2;
  auto ff = [&f, scale, mean](double x) { return f(scale * x + mean); };

  const double r1 = gauss_kronrod_21_m1_1(ff, &error_local, l1);
  double estimate = scale * r1;

  const double abs_tol1 = std::abs(estimate * relative_tolerance);
  if (is_root) {
    absolute_tolerance = std::max(absolute_tolerance, abs_tol1);
  }

  if (max_levels && (abs_tol1 < error_local)
      && (absolute_tolerance < error_local)) {
    const double mid = (a + b) / 2;
    double l1_local = 0.0;
    estimate = gauss_kronrod_21_recursive(
        f, relative_tolerance, a, mid, max_levels - 1, absolute_tolerance / 2,
        false, error, l1);
    estimate += gauss_kronrod_21_recursive(
        f, relative_tolerance, mid, b, max_levels - 1, absolute_tolerance / 2,
        false, &error_local, &l1_local);
    *error += error_local;
    *l1 += l1_local;
    return estimate;
  }

  *l1 *= scale;
  *error = error_local;
  return estimate;
}

/**
 * Integrate f from a to b with adaptive (G10, K21) quadrature, honouring an
 * absolute-error floor during refinement.
 *
 * Finite, reversed, semi-infinite and doubly-infinite limits use Boost's
 * changes of variable.  The absolute floor is passed through those changes of
 * variable unscaled, which keeps it in the same units as the error written to
 * *error (Boost leaves both in reference-panel units; see
 * gauss_kronrod_21_m1_1).
 *
 * With absolute_tolerance == 0 this is bit-for-bit equivalent to
 * boost::math::quadrature::gauss_kronrod<double, 21>::integrate with the same
 * arguments.  The prim test
 * StanMath_integrate_1d_gk_prim.matches_boost_bit_for_bit_when_abs_tol_zero
 * pins that equivalence.
 *
 * @tparam F type of f
 * @param f function to integrate
 * @param a lower limit of integration (may be -infinity)
 * @param b upper limit of integration (may be +infinity)
 * @param max_depth maximum bisection depth
 * @param relative_tolerance local relative target per panel
 * @param absolute_tolerance absolute error floor on refinement
 * @param[out] error estimate of the quadrature error (must not be null)
 * @param[out] l1 estimate of the L1 norm (must not be null)
 * @return estimate of the integral of f from a to b
 * @throw std::domain_error if the limits are NaN or otherwise not sensible
 */
template <typename F>
inline double gauss_kronrod_21_integrate(const F& f, double a, double b,
                                         unsigned int max_depth,
                                         double relative_tolerance,
                                         double absolute_tolerance,
                                         double* error, double* l1) {
  static constexpr const char* function = "gauss_kronrod_21_integrate";
  const double max_value = boost::math::tools::max_value<double>();

  if (!std::isnan(a) && !std::isnan(b)) {
    if (a <= -max_value && b >= max_value) {
      auto u = [&f](double t) {
        const double t_sq = t * t;
        const double inv = 1 / (1 - t_sq);
        const double w = (1 + t_sq) * inv * inv;
        return f(t * inv) * w;
      };
      return gauss_kronrod_21_recursive(u, relative_tolerance, -1.0, 1.0,
                                        max_depth, absolute_tolerance, true,
                                        error, l1);
    }

    if (std::isfinite(a) && b >= max_value) {
      auto u = [&f, a](double t) {
        const double z = 1 / (t + 1);
        return f(2 * z + a - 1) * z * z;
      };
      const double Q = 2
                       * gauss_kronrod_21_recursive(
                           u, relative_tolerance, -1.0, 1.0, max_depth,
                           absolute_tolerance, true, error, l1);
      *l1 *= 2;
      return Q;
    }

    if (std::isfinite(b) && a <= -max_value) {
      auto v = [&f, b](double t) {
        const double z = 1 / (t + 1);
        return f(b - (2 * z - 1)) * z * z;
      };
      const double Q = 2
                       * gauss_kronrod_21_recursive(
                           v, relative_tolerance, -1.0, 1.0, max_depth,
                           absolute_tolerance, true, error, l1);
      *l1 *= 2;
      return Q;
    }

    if (std::isfinite(a) && std::isfinite(b)) {
      if (a == b) {
        *error = 0.0;
        *l1 = 0.0;
        return 0.0;
      }
      if (b < a) {
        return -gauss_kronrod_21_recursive(f, relative_tolerance, b, a,
                                           max_depth, absolute_tolerance, true,
                                           error, l1);
      }
      return gauss_kronrod_21_recursive(f, relative_tolerance, a, b, max_depth,
                                        absolute_tolerance, true, error, l1);
    }
  }

  throw_domain_error(function, "integration limits", a, "lower limit is ",
                     ", and the domain of integration is not sensible");
  return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
