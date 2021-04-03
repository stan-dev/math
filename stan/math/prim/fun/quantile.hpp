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
  check_nonzero_size("quantile", "xs", xs);
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

template <typename T, typename Tp, require_all_vector_t<T, Tp>* = nullptr>
inline Eigen::ArrayXd quantile(const T& xs, const Tp& ps) {
  check_not_nan("quantile", "container argument", xs);
  check_bounded("quantile", "ps", ps, 0, 1);
  check_nonzero_size("quantile", "xs", xs);
  check_nonzero_size("quantile", "ps", ps);
  size_t n_sample = xs.size();
  size_t n_ps = ps.size();

  Eigen::VectorXd x = as_array_or_scalar(xs);
  Eigen::ArrayXd p = as_array_or_scalar(ps);
  Eigen::ArrayXd ret = Eigen::ArrayXd::Zero(n_ps);

  if (n_sample == 1)
    return ret + x[0];

  std::sort(x.data(), x.data() + n_sample);
  Eigen::ArrayXd index = (n_sample - 1) * p;

  for (size_t i = 0; i < n_ps; ++i) {
    if (p[i] == 0.)
      ret[i] = *std::min_element(x.data(), x.data() + n_sample);
    if (p[i] == 1.)
      ret[i] = *std::max_element(x.data(), x.data() + n_sample);

    size_t lo = std::floor(index[i]);
    size_t hi = std::ceil(index[i]);

    double h = index[i] - lo;

    ret[i] = (1 - h) * x.coeff(lo) + h * x.coeff(hi);
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
