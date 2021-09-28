#ifndef STAN_MATH_PRIM_PROB_DISCRETE_RANGE_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_DISCRETE_RANGE_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return the log CCDF of a discrete range distribution for the given y,
 * lower and upper bounds (all integers).
 *
 * `y`, `lower` and `upper` can each be a scalar or a one-dimensional container.
 * Any container arguments must be the same size.
 *
 * @tparam T_y type of scalar, either int or std::vector<int>
 * @tparam T_lower type of lower bound, either int or std::vector<int>
 * @tparam T_upper type of upper bound, either int or std::vector<int>
 *
 * @param y integer random variable
 * @param lower integer lower bound
 * @param upper integer upper bound
 * @return The log CCDF evaluated at the specified arguments. If containers are
 * supplied, returns the sum of the log CCDFs.
 * @throw std::domain_error if upper is smaller than lower.
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes.
 */
template <typename T_y, typename T_lower, typename T_upper>
double discrete_range_lccdf(const T_y& y, const T_lower& lower,
                            const T_upper& upper) {
  static const char* function = "discrete_range_lccdf";
  check_consistent_sizes(function, "Lower bound parameter", lower,
                         "Upper bound parameter", upper);
  check_greater_or_equal(function, "Upper bound parameter", upper, lower);

  if (size_zero(y, lower, upper)) {
    return 0;
  }

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_lower> lower_vec(lower);
  scalar_seq_view<T_upper> upper_vec(upper);
  size_t N = max_size(y, lower, upper);

  for (size_t n = 0; n < N; ++n) {
    const int y_dbl = y_vec[n];
    if (y_dbl >= upper_vec[n]) {
      return LOG_ZERO;
    }
  }

  double ccdf(0.0);
  for (size_t n = 0; n < N; n++) {
    const int y_dbl = y_vec[n];
    if (y_dbl >= lower_vec[n]) {
      const int lower_dbl = lower_vec[n];
      const int upper_dbl = upper_vec[n];
      ccdf += log(upper_dbl - y_dbl) - log(upper_dbl - lower_dbl + 1);
    }
  }
  return ccdf;
}

}  // namespace math
}  // namespace stan
#endif
