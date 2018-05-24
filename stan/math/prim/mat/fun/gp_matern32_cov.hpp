#ifndef STAN_MATH_PRIM_MAT_FUN_GP_MATERN32_COV_HPP
#define STAN_MATH_PRIM_MAT_FUN_GP_MATERN32_COV_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a Matern 3/2 covariance matrix with one input vector
 *
 * \f[ k(x, x') = \sigma^2(1 +
 *  \frac{\sqrt{3}d(x, x')}{l})exp(-\frac{\sqrt{3}d(x, x')}{l})
 * \f]
 *
 * where \f$ d(x, x') \f$ is the Euclidean distance.
 *
 * @tparam T_x type of elements contained in vector x
 * @tparam T_s type of element of sigma, marginal standard deviation
 * @tparam T_l type of elements of length scale
 * 
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern32_cov(const std::vector<T_x> &x, const T_s &sigma,
                  const T_l &length_scale) {
  using std::exp;
  using std::pow;

  size_t x_size = x.size();
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern32_cov", "x", x[n]);

  check_positive_finite("gp_matern32_cov", "marginal variance", sigma);
  check_positive_finite("gp_matern32_cov", "length scale", length_scale);

  Eigen::Matrix<typename return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3_inv_l = sqrt(3.0) / length_scale;
  T_l neg_root_3_inv_l = -1.0 * sqrt(3.0) / length_scale;
  typename return_type<T_x>::type distance;
  
  for (size_t i = 0; i < x_size; ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      distance = sqrt(squared_distance(x[i], x[j]));
      cov(i, j) = sigma_sq * (1.0 + root_3_inv_l * distance) *
        exp(neg_root_3_inv_l * distance);
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}

/**
 * Returns a Matern 3/2 Kernel with one input vector,
 * with automatic relevance determination (ARD)
 *
 * \f[ k(x, x') = \sigma^2(1 + \sqrt{3}
 *   \sqrt{\sum_{k=1}^{K}\frac{d{(x, x')^2}}{l_k^2}})
 *   exp(-\sqrt{3}\sqrt{\sum_{k=1}^{K}\frac{d(x, x')^2}{l_k^2}}) \f]
 *
 * where \f$d(x, x')\f$ is the Euclidean distance.
 *
 * @tparam T_x type of elements contained in vector x
 * @tparam T_s type of element of sigma, marginal standard deviation
 * @tparam T_l type of elements of length scale
 * 
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern32_cov(const std::vector<T_x> &x, const T_s &sigma,
                  const std::vector<T_l> &length_scale) {
  using std::exp;

  size_t x_size = x.size();
  size_t l_size = length_scale.size();
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern32_cov", "x", x[n]);

  check_positive_finite("gp_matern32_cov", "marginal variance", sigma);
  check_positive_finite("gp_matern32_cov", "length scale", length_scale);

  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern32_cov", "length scale", length_scale[n]);

  Eigen::Matrix<typename return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3 = sqrt(3.0);
  T_l neg_root_3 = -1.0 * sqrt(3.0);
  T_l temp;
  typename return_type<T_x>::type sq_distance;
  
  for (size_t i = 0; i < x_size; ++i) {
    for (size_t j = i; j < x_size; ++j) {
      temp = 0;
      sq_distance = squared_distance(x[i], x[j]);
      for (size_t k = 0; k < l_size; ++k) {
        temp += 1.0 / square(length_scale[k]);
      }
      cov(i, j) = sigma_sq * (1.0 + root_3 * sqrt(sq_distance * temp)) *
        exp(neg_root_3 * sqrt(sq_distance * temp));
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}

/**
 * Returns a Matern 3/2 covariance matrix with two input vectors
 *
 * \f[ k(x, x') = \sigma^2(1 +
 *  \frac{\sqrt{3}d(x, x')}{l})exp(-\sqrt{3}\frac{d(x, x')}{l})
 * \f]
 *
 * where \f$d(x, x')\f$ is the Euclidean distance.
 *
 * @tparam T_x1 type of elements contained in vector x1
 * @tparam T_x2 type of elements contained in vector x2
 * @tparam T_s type of element of sigma, marginal standard deviation
 * @tparam T_l type of elements of length scale
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
    typename return_type<T_x1, T_x2, T_s, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_matern32_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                  const T_s &sigma, const T_l &length_scale) {
  using std::exp;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern32_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern32_cov", "x2", x2[n]);

  check_positive_finite("gp_matern32_cov", "marginal variance", sigma);
  check_positive_finite("gp_matern32_cov", "length scale", length_scale);

  Eigen::Matrix<typename return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3_inv_l_sq = sqrt(3.0) / length_scale;
  T_l neg_root_3_inv_l_sq = -1.0 * sqrt(3.0) / length_scale;
  typename return_type<T_x1, T_x2>::type distance;
  
  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      distance = sqrt(squared_distance(x1[i], x2[j]));
      cov(i, j) = sigma_sq * (1.0 + root_3_inv_l_sq * distance) *
        exp(neg_root_3_inv_l_sq * distance);
    }
  }
  return cov;
}

/**
 * Returns a Matern 3/2 Kernel with two input vectors with automatic
 * relevance determination (ARD)
 *
 * \f[ k(x, x') = \sigma^2(1 + \sqrt{3}
 *   \sqrt{\sum_{k=1}^{K}\frac{d(x, x')^}{l_k^2}})
 *   exp(-\sqrt{3}\sqrt{\sum_{k=1}^{K}\frac{d(x, x')^2}{l_k^2}}) \f]
 *
 * where \f$d(x, x')\f$ is the Euclidean distance
 *
 * @tparam T_x1 type of elements contained in vector x1
 * @tparam T_x2 type of elements contained in vector x2
 * @tparam T_s type of element of sigma, marginal standard deviation
 * @tparam T_l type of elements of length scale
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
    typename return_type<T_x1, T_x2, T_s, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_matern32_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                  const T_s &sigma, const std::vector<T_l> &length_scale) {
  using std::exp;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  size_t l_size = length_scale.size();

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern32_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern32_cov", "x2", x2[n]);

  check_positive_finite("gp_matern32_cov", "marginal variance", sigma);
  check_positive_finite("gp_matern32_cov", "length scale", length_scale);

  for (size_t n = 0; n < l_size; ++n)
    check_not_nan("gp_matern32_cov", "length scale", length_scale[n]);

  Eigen::Matrix<typename return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3 = sqrt(3.0);
  T_l neg_root_3 = -1.0 * sqrt(3.0);
  T_l temp;
  typename return_type<T_x1, T_x2>::type sq_distance;
  
  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      temp = 0;
      sq_distance = squared_distance(x1[i], x2[j]);
      for (size_t k = 0; k < l_size; ++k) {
        temp += 1.0 / square(length_scale[k]);
      }
      cov(i, j) = sigma_sq * (1.0 + root_3 * sqrt(sq_distance * temp)) *
        exp(neg_root_3 * sqrt(sq_distance * temp));
    }
  }
  return cov;
}
}  // namespace math
}  // namespace stan
#endif
