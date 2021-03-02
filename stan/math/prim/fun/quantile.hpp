#ifndef STAN_MATH_PRIM_FUN_QUANTILE_HPP
#define STAN_MATH_PRIM_FUN_QUANTILE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * Return sample quantiles corresponding to the given probabilities.
 * The smallest observation corresponds to a probability of 0 and the largest to
 * a probability of 1.
 *
 * Implementation follows the default R behavior, as defined here:
 * https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/quantile
 *
 * @tparam T Type of elements contained in vector.
 * @param xs Numeric vector whose sample quantiles are wanted
 * @param p Probability with value between 0 and 1.
 * @return Sample quantile.
 * @throw std::domain_error If any of the values are NaN.
 */
template <typename T, require_vector_t<T>* = nullptr>
inline double quantile(const T& xs, const double p) {
  using boost::accumulators::accumulator_set;
  using boost::accumulators::left;
  using boost::accumulators::quantile;
  using boost::accumulators::quantile_probability;
  using boost::accumulators::right;
  using boost::accumulators::stats;
  using boost::accumulators::tag::tail;
  using boost::accumulators::tag::tail_quantile;

  check_not_nan("quantile", "container argument", xs);
  check_bounded("quantile", "p", p, 0, 1);

  const auto& x = as_array_or_scalar(xs);

  double M = x.rows();

  size_t cache_size = M;

  if (p < 0.5) {
    accumulator_set<double, stats<tail_quantile<left> > > acc(
        tail<left>::cache_size = cache_size);
    for (int i = 0; i < M; i++)
      acc(x(i));
    return quantile(acc, quantile_probability = p);
  }
  accumulator_set<double, stats<tail_quantile<right> > > acc(
      tail<right>::cache_size = cache_size);
  for (int i = 0; i < M; i++)
    acc(x(i));
  return quantile(acc, quantile_probability = p);
}

}  // namespace math
}  // namespace stan

#endif
