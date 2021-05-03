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
 * @tparam T An Eigen type with either fixed rows or columns at compile time, or
 * std::vector<double>
 * @param samples_vec Numeric vector whose sample quantiles are wanted
 * @param p Probability with value between 0 and 1.
 * @return Sample quantile.
 * @throw std::invalid_argument If any element of samples_vec is NaN, or size 0.
 * @throw std::domain_error If `p<0` or `p>1`.
 */
template <typename T, require_vector_t<T>* = nullptr,
          require_vector_vt<std::is_arithmetic, T>* = nullptr>
inline double quantile(const T& samples_vec, const double p) {
  check_not_nan("quantile", "p", p);
  check_bounded("quantile", "p", p, 0, 1);

  const size_t n_sample = samples_vec.size();
  if (n_sample == 0) {
    return {};
  }

  Eigen::VectorXd x = as_array_or_scalar(samples_vec);
  check_not_nan("quantile", "samples_vec", x);

  if (n_sample == 1)
    return x.coeff(0);
  else if (p == 0.)
    return *std::min_element(x.data(), x.data() + n_sample);
  else if (p == 1.)
    return *std::max_element(x.data(), x.data() + n_sample);

  double index = (n_sample - 1) * p;
  size_t lo = std::floor(index);
  size_t hi = std::ceil(index);

  std::sort(x.data(), x.data() + n_sample, std::less<double>());

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
 * @tparam T An Eigen type with either fixed rows or columns at compile time, or
 * std::vector<double>
 * @tparam Tp Type of probabilities vector
 * @param samples_vec Numeric vector whose sample quantiles are wanted
 * @param ps Vector of probability with value between 0 and 1.
 * @return Sample quantiles, one for each p in ps.
 * @throw std::invalid_argument If any of the values are NaN or size 0.
 * @throw std::domain_error If `p<0` or `p>1` for any p in ps.
 */
template <typename T, typename Tp, require_all_vector_t<T, Tp>* = nullptr,
          require_vector_vt<std::is_arithmetic, T>* = nullptr,
          require_std_vector_vt<std::is_arithmetic, Tp>* = nullptr>
inline std::vector<double> quantile(const T& samples_vec, const Tp& ps) {
  check_not_nan("quantile", "ps", ps);
  check_bounded("quantile", "ps", ps, 0, 1);

  const size_t n_sample = samples_vec.size();
  const size_t n_ps = ps.size();
  if (n_ps == 0 || n_sample == 0) {
    return {};
  }

  Eigen::VectorXd x = as_array_or_scalar(samples_vec);
  check_not_nan("quantile", "samples_vec", x);

  const auto& p = as_array_or_scalar(ps);
  std::vector<double> ret(n_ps, 0.0);

  std::sort(x.data(), x.data() + n_sample, std::less<double>());
  Eigen::ArrayXd index = (n_sample - 1) * p;

  for (size_t i = 0; i < n_ps; ++i) {
    if (p[i] == 0.) {
      ret[i] = x.coeff(0);
    } else if (p[i] == 1.) {
      ret[i] = x.coeff(n_sample - 1);
    } else {
      size_t lo = std::floor(index[i]);
      size_t hi = std::ceil(index[i]);

      double h = index[i] - lo;

      ret[i] = (1 - h) * x.coeff(lo) + h * x.coeff(hi);
    }
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
