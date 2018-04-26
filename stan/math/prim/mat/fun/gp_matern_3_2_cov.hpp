#ifndef STAN_MATH_PRIM_MAT_FUN_GP_MATERN_3_2_COV_HPP
#define STAN_MATH_PRIM_MAT_FUN_GP_MATERN_3_2_COV_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/** Returns a Matern 3/2 covariance matrix with one input vector
 *
 * \f[ k(x, x') = \sigma^2(1 + \sqrt{3}
 *  \frac{\sqrt{(x - x')^2}}{l})exp(-\sqrt{3}\frac{\sqrt{(x - x')^2}}{l})
 * \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern_3_2_cov(const std::vector<T_x> &x, const T_s &sigma,
                  const T_l &length_scale) {
  using stan::math::squared_distance;
  using std::abs;
  using std::exp;
  using std::pow;

  size_t x_size = x.size();

  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);

  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x", x[n]);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3_inv_l = pow(3, 0.5) / length_scale;
  T_l neg_root_3_inv_l = -1.0 * pow(3, 0.5) / length_scale;

  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_sq *
                  (1.0 + root_3_inv_l * squared_distance(x[i], x[j])) *
                  exp(neg_root_3_inv_l * squared_distance(x[i], x[j]));
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}

/** Returns a Matern 3/2 Kernel with one input vector,
 * with automatic relevance determination (ARD) priors
 *
 * \f[ k(x, x') = \sigma^2(1 + \sqrt{3}
 *   \frac{\sum_{k=1}^{K}\sqrt{(x - x')^2}}{l_k}
 *   exp(-\sqrt{3}\sum_{k=1}^{K}\frac{\sqrt{(x - x')^2}}{l_k}) \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern_3_2_cov(const std::vector<T_x> &x, const T_s &sigma,
                  const std::vector<T_l> &length_scale) {
  using std::abs;
  using std::exp;
  using std::pow;

  size_t x_size = x.size();
  size_t l_size = length_scale.size();

  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);

  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x", x[n]);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3 = pow(3.0, 0.5);
  T_l neg_root_3 = -1.0 * pow(3.0, 0.5);
  T_l temp;

  for (size_t i = 0; i < x_size; ++i) {
    for (size_t j = i; j < x_size; ++j) {
      temp = 0;
      for (size_t k = 0; k < l_size; ++k) {
        temp += squared_distance(x[i], x[j]) / length_scale[k];
      }
      temp = pow(temp, 0.5);
      cov(i, j) = sigma_sq * (1.0 + root_3 * temp) * exp(neg_root_3 * temp);
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}

/** Returns a Matern 3/2 covariance matrix with two input vectors
 *
 * \f[ k(x, x') = \sigma^2(1 + \sqrt{3}
 *  \frac{\sqrt{(x - x')^2}}{l}exp(-\sqrt{3}\frac{\sqrt{(x - x')^2}}{l})
 * \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 *
 * @param x1 std::vector of elements that can be used in stan::math::distance
 * @param x2 std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x1, x2 are nan or inf
 *
 */
template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_s, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_matern_3_2_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                  const T_s &sigma, const T_l &length_scale) {
  using stan::math::squared_distance;
  using std::abs;
  using std::exp;
  using std::pow;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();

  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x2", x2[n]);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3_inv_l = pow(3.0, 0.5) / length_scale;
  T_l neg_root_3_inv_l = -1.0 * pow(3.0, 0.5) / length_scale;

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = sigma_sq *
                  (1.0 + root_3_inv_l * squared_distance(x1[i], x2[j])) *
                  exp(neg_root_3_inv_l * squared_distance(x1[i], x2[j]));
    }
  }
  return cov;
}

/** Returns a Matern 3/2 Kernel with two input vectors with automatic
 * relevance determination (ARD) priors
 *
 * \f[ k(x, x') = \sigma^2(1 + \sqrt{3}
 *   \frac{\sum_{k=1}^{K}\sqrt{(x - x')^2}}{l_k}
 *   exp(-\sqrt{3}\sum_{k=1}^{K}\frac{\sqrt{(x - x')^2}}{l_k}) \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 *
 * @param x1 std::vector of elements that can be used in stan::math::distance
 * @param x2 std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x1, x2 are nan or inf
 *
 */
template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_s, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_matern_3_2_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                  const T_s &sigma, const std::vector<T_l> &length_scale) {
  using stan::math::squared_distance;
  using std::abs;
  using std::exp;
  using std::pow;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  size_t l_size = length_scale.size();

  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x2", x2[n]);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3 = pow(3.0, 0.5);
  T_l neg_root_3 = -1.0 * pow(3.0, 0.5);
  T_l temp;

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      temp = 0;
      for (size_t k = 0; k < l_size; ++k) {
        temp += squared_distance(x1[i], x2[j]) / length_scale[k];
      }
      temp = pow(temp, 0.5);
      cov(i, j) = sigma_sq * (1.0 + root_3 * temp) * exp(neg_root_3 * temp);
    }
  }
  return cov;
}
} // namespace math
} // namespace stan
#endif
