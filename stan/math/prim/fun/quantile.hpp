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
 * Implements algorithm 7 from Hyndman, R. J. and Fan, Y.,
 * Sample quantiles in Statistical Packages (R's default quantile function)
 *
 * @tparam T Type of input vector
 * @param xs Numeric vector whose sample quantiles are wanted
 * @param p Probability with value between 0 and 1.
 * @return Sample quantile.
 * @throw std::invalid_argument If any element of xs is NaN, or size 0.
 * @throw std::domain_error If p<0 or p>1.
 */
template <typename T, require_vector_t<T>* = nullptr,
          require_vector_vt<std::is_arithmetic, T>* = nullptr>
inline double quantile(const T& xs, const double p) {
  check_not_nan("quantile", "xs", xs);
  check_nonzero_size("quantile", "xs", xs);

  check_not_nan("quantile", "p", p);
  check_bounded("quantile", "p", p, 0, 1);

  size_t n_sample = xs.size();

  Eigen::VectorXd x = as_array_or_scalar(xs);

  if (n_sample == 1)
    return x.coeff(0);
  if (p == 0.)
    return *std::min_element(x.data(), x.data() + n_sample);
  if (p == 1.)
    return *std::max_element(x.data(), x.data() + n_sample);

  double index = (n_sample - 1) * p;
  size_t lo = std::floor(index);
  size_t hi = std::ceil(index);

  std::sort(x.data(), x.data() + n_sample);

  double h = index - lo;
  return (1 - h) * x.coeff(lo) + h * x.coeff(hi);
}

/**
 * Return sample quantiles corresponding to the given probabilities.
 * The smallest observation corresponds to a probability of 0 and the largest to
 * a probability of 1.
 *
 * Implements algorithm 7 from Hyndman, R. J. and Fan, Y.,
 * Sample quantiles in Statistical Packages (R's default quantile function)
 *
 * @tparam T Type of input vector
 * @tparam Tp Type of probabilities vector
 * @param xs Numeric vector whose sample quantiles are wanted
 * @param ps Vector of probability with value between 0 and 1.
 * @return Sample quantiles, one for each p in ps.
 * @throw std::invalid_argument If any of the values are NaN or size 0.
 * @throw std::domain_error If p<0 or p>1 for any p in ps.
 */
template <typename T, typename Tp, require_all_vector_t<T, Tp>* = nullptr,
          require_vector_vt<std::is_arithmetic, T>* = nullptr,
          require_vector_vt<std::is_arithmetic, Tp>* = nullptr>
inline std::vector<double> quantile(const T& xs, const Tp& ps) {
  check_not_nan("quantile", "xs", xs);
  check_not_nan("quantile", "ps", ps);

  check_nonzero_size("quantile", "xs", xs);
  check_nonzero_size("quantile", "ps", ps);

  check_bounded("quantile", "ps", ps, 0, 1);

  size_t n_sample = xs.size();
  size_t n_ps = ps.size();

  Eigen::VectorXd x = as_array_or_scalar(xs);
  const auto& p = as_array_or_scalar(ps);
  std::vector<double> ret(n_ps, 0.0);

  std::sort(x.data(), x.data() + n_sample);
  Eigen::ArrayXd index = (n_sample - 1) * p;

  for (size_t i = 0; i < n_ps; ++i) {
    if (p[i] == 0.)
      ret[i] = x.coeff(0);
    if (p[i] == 1.)
      ret[i] = x.coeff(n_sample - 1);

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
