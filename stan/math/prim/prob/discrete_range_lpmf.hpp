#ifndef STAN_MATH_PRIM_PROB_DISCRETE_RANGE_LPMF_HPP
#define STAN_MATH_PRIM_PROB_DISCRETE_RANGE_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return the log PMF of a discrete range for the given y, lower and upper
 * bound (all integers).
 *
 \f{eqnarray*}{
   y &\sim& \mbox{\sf{discrete\_range}}(lower, upper) \\
     \log(p (y \, |\, lower, upper))
        &=& \log \left( \frac{1}{upper - lower + 1} \right) \\
        &=& -\log (upper - lower + 1)
 \f}
 *
 * `lower` and `upper` can each be a scalar or a one-dimensional container.
 * Any container arguments must be the same size.
 *
 * @tparam T_y type of scalar, either int or std::vector<int>
 * @tparam T_lower type of lower bound, either int or std::vector<int>
 * @tparam T_upper type of upper bound, either int or std::vector<int>
 *
 * @param y integer random variable
 * @param lower integer lower bound
 * @param upper integer upper bound
 * @return Log probability. If containers are supplied, returns the log sum
 * of the probabilities.
 * @throw std::domain_error if upper is smaller than lower.
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes.
 */
template <bool propto, typename T_y, typename T_lower, typename T_upper>
double discrete_range_lpmf(const T_y& y, const T_lower& lower,
                           const T_upper& upper) {
  using std::log;
  static constexpr const char* function = "discrete_range_lpmf";
  check_not_nan(function, "Random variable", y);
  check_consistent_sizes(function, "Lower bound parameter", lower,
                         "Upper bound parameter", upper);
  check_greater_or_equal(function, "Upper bound parameter", upper, lower);

  if (size_zero(y, lower, upper)) {
    return 0.0;
  }
  if (!include_summand<propto>::value) {
    return 0.0;
  }

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_lower> lower_vec(lower);
  scalar_seq_view<T_upper> upper_vec(upper);
  size_t size_lower_upper = max_size(lower, upper);
  size_t N = max_size(y, lower, upper);

  for (size_t n = 0; n < N; ++n) {
    const double y_dbl = y_vec.val(n);
    if (y_dbl < lower_vec.val(n) || y_dbl > upper_vec.val(n)) {
      return LOG_ZERO;
    }
  }

  VectorBuilder<true, double, T_lower, T_upper> log_upper_minus_lower(
      size_lower_upper);

  for (size_t i = 0; i < size_lower_upper; i++) {
    const double lower_dbl = lower_vec.val(i);
    const double upper_dbl = upper_vec.val(i);
    log_upper_minus_lower[i] = log(upper_dbl - lower_dbl + 1);
  }

  double logp(0.0);
  for (size_t n = 0; n < N; n++) {
    logp -= log_upper_minus_lower[n];
  }

  return logp;
}

template <typename T_y, typename T_lower, typename T_upper>
inline double discrete_range_lpmf(const T_y& y, const T_lower& lower,
                                  const T_upper& upper) {
  return discrete_range_lpmf<false>(y, lower, upper);
}

}  // namespace math
}  // namespace stan
#endif
