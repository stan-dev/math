#ifndef STAN_MATH_PRIM_FUN_QUANTILE_HPP
#define STAN_MATH_PRIM_FUN_QUANTILE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
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
  check_not_nan("quantile", "container argument", xs);
  check_bounded("quantile", "p", p, 0, 1);
  size_t n_sample = xs.size();

  Eigen::VectorXd x = as_array_or_scalar(xs);

  if (n_sample == 1)
    return x[0];
  if (p == 0.)
    return *std::min_element(x.data(), x.data() + n_sample);
  if (p == 1.)
    return *std::max_element(x.data(), x.data() + n_sample);

  double index = (n_sample - 1) * p;
  size_t lo = std::floor(index);
  size_t hi = std::ceil(index);

  std::sort(x.data(), x.data() + n_sample);

  double q = x[lo];
  double h = index - lo;
  return (1 - h) * q + h * x[hi];
}

}  // namespace math
}  // namespace stan

#endif
