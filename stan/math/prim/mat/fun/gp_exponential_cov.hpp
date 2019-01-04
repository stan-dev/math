#ifndef STAN_MATH_PRIM_MAT_FUN_GP_EXPONENTIAL_COV_HPP
#define STAN_MATH_PRIM_MAT_FUN_GP_EXPONENTIAL_COV_HPP

#include <cmath>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/distance.hpp>
#include <stan/math/prim/mat/fun/divide_columns.hpp>
#include <stan/math/prim/mat/fun/sqrt.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a Matern exponential covariance matrix with one input vector
 *
 * \f[ k(x, x') = \sigma^2 exp(-\frac{d(x, x')^2}{2l^2}) \f]
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param sigma standard deviation that can be used in stan::math::square
 * @param length_scale length scale
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_exponential_cov(const std::vector<T_x> &x, const T_s &sigma,
                   const T_l &length_scale) {
  using std::exp;
  using std::pow;

  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i)
    check_not_nan("gp_exponential_cov", "x", x[i]);

  check_positive_finite("gp_exponential_cov", "magnitude", sigma);
  check_positive_finite("gp_exponential_cov", "length scale", length_scale);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);
  
  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l neg_inv_l = -1.0 / length_scale;

  for (size_t i = 0; i < x_size; ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_sq * exp(neg_inv_l * squared_distance(x[i], x[j]));
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}

/**
 * Returns a Matern exponential covariance matrix with one input vector
 *  with automatic relevance determination (ARD)
 *
 * \f[ k(x, x') = \sigma^2 exp(-\sum_{k=1}^K\frac{d(x, x')}{l_k}) \f]
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param sigma standard deviation that can be used in stan::math::square
 * @param length_scale length scale
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_exponential_cov(const std::vector<Eigen::Matrix<T_x, -1, 1>> &x,
                   const T_s &sigma, const std::vector<T_l> &length_scale) {
  using std::exp;
  using std::pow;

  size_t x_size = x.size();
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_exponential_cov", "x", x[n]);

  check_positive_finite("gp_exponential_cov", "marginal variance", sigma);
  check_positive_finite("gp_exponential_cov", "length scale", length_scale);

  size_t l_size = length_scale.size();
  check_size_match("gp_exponential_cov", "x dimension", x[0].size(),
                   "number of length scales", l_size);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  std::vector<T_x> x_new = divide_columns(x, length_scale);

  T_s sigma_sq = square(sigma);
  double neg_one_half = -1.0 / 2.0;
  // for (size_t i = 0; i < (x_size - 1); ++i) {
  //   cov(i, i) = sigma_sq;
  //   for (size_t j = i + 1; j < x_size; ++j) {
  //     typename return_type<T_x, T_l>::type sq_distance =
  //       squared_distance(x_new[i], x_new[j]);
  //     cov(i, j) =
  //         sigma_sq * exp(neg_one_half * sq_distance);
  //     cov(j, i) = cov(i, j);
  //   }
  // }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}


 * Returns a Matern exponential covariance matrix
 *
 * \f[ k(x, x') = \sigma^2 exp(-\sum_{k=1}^K\frac{d(x, x')}{l_k}) \f]
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param sigma standard deviation that can be used in stan::math::square
 * @param length_scale length scale
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_exponential_cov(const std::vector<Eigen::Matrix<T_x, -1, 1>> &x,
                   const T_s &sigma, const T_l &length_scale) {
  using std::exp;
  using std::pow;

  size_t x_size = x.size();
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_exponential_cov", "x", x[n]);

  check_positive_finite("gp_exponential_cov", "marginal variance", sigma);
  check_positive_finite("gp_exponential_cov", "length scale", length_scale);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  double neg_one_half = -1.0 / 2.0;
  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      typename return_type<T_x, T_l>::type sq_distance =
        squared_distance(x_new[i], x_new[j]);
      cov(i, j) =
          sigma_sq * exp(neg_one_half * sq_distance);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}
  
/**
 * Returns a Matern exponential covariance matrix with two input vectors
 *
 * \f[ k(x, x') = \sigma^2 exp(-\frac{d(x, x')}{l}) \f]
 *
 * @param x1 std::vector of elements that can be used in stan::math::distance
 * @param x2 std::vector of elements that can be used in stan::math::distance
 * @param sigma standard deviation that can be used in stan::math::square
 * @param length_scale length scale
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_s, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_exponential_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                   const T_s &sigma, const T_l &length_scale) {
  using std::exp;
  using std::pow;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_exponential_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_exponential_cov", "x2", x2[n]);

  check_positive_finite("gp_exponential_cov", "marginal variance", sigma);
  check_positive_finite("gp_exponential_cov", "length scale", length_scale);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l neg_inv_l = -1.0 / length_scale;

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) =
          sigma_sq * exp(neg_inv_l * squared_distance(x1[i], x2[j]));
    }
  }
  return cov;
}

/**
 * Returns a Matern exponential covariance matrix with two input vectors
 *  with automatic relevance determination (ARD)
 *
 * \f[ k(x, x') = \sigma^2 exp(-\sum_{k=1}^K\frac{d(x, x')}{l_k}) \f]
 *
 * @param x1 std::vector of elements that can be used in stan::math::distance
 * @param x2 std::vector of elements that can be used in stan::math::distance
 * @param sigma standard deviation that can be used in stan::math::square
 * @param length_scale length scale
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_s, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_exponential_cov(const std::vector<Eigen::Matrix<T_x1, -1, 1>> &x1,
                   const std::vector<Eigen::Matrix<T_x2, -1, 1>> &x2,
                   const T_s &sigma, const std::vector<T_l> &length_scale) {
  using std::exp;
  using std::pow;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  size_t l_size = length_scale.size();

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_exponential_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_exponential_cov", "x2", x2[n]);

  check_positive_finite("gp_exponential_cov", "marginal variance", sigma);
  check_positive_finite("gp_exponential_cov", "length scale", length_scale);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l temp;

  // fix ARD part

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) =
          sigma_sq * exp(-1.0 * sqrt(squared_distance(x1[i], x2[j])) * temp);
    }
  }
  return cov;
}
}  // namespace math
}  // namespace stan
#endif
