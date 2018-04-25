#ifndef STAN_MATH_PRIM_MAT_FUN_GP_MATERN_3_2_COV_HPP
#define STAN_MATH_PRIM_MAT_FUN_GP_MATERN_3_2_COV_HPP

#include <cmath>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/fun/squared_distance.hpp>

namespace stan {
namespace math {

/** Returns a Matern 3/2 covariance matrix with one input vector
 *
 * \f[ k(x, x') = \sigma^2(1 + \gamma \sqrt{3}
 *  \frac{\sqrt{(x - x')^2}}{l}exp(-\gamma\sqrt{3}\frac{\sqrt{(x - x')^2}}{l})  \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 * 
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @param gamma
 * @throw std::domain error if sigma <= 0, l <= 0, or x or gamma are nan of inf
 *
 */
template <typename T_x, typename T_l, typename T_s, typename T_g>
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_l, T_s, T_g>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    gp_matern_3_2_cov(const std::vector<T_x> &x, const T_l &length_scale,
                      const T_s &sigma = 1.0, const T_g &gamma = 1.0) {
  using std::pow;
  using std::abs;
  using std::exp;
  using stan::math::squared_distance;

  size_t x_size = x.size();
  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);

  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_not_nan("gp_matern_3_2_cov", "gamma", gamma);
  
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x", x[n]);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_g, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3_inv_l = pow(3, 0.5) / length_scale;
  T_g neg_gamma = -1.0 * gamma;

  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_sq * (1.0 + root_3_inv_l * gamma *
                              squared_distance(x[i], x[j])) *
        exp(neg_gamma * root_3_inv_l * squared_distance(x[i], x[j]));
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}

/** Returns a Matern 3/2 Kernel with one input vector,
 * with automatic relevance determination (ARD) priors
 *
 * \f[ k(x, x') = \sigma^2(1 + \gamma \sqrt{3} 
 *   \frac{\sum_{k=1}^{K}\sqrt{(x - x')^2}}{l_k}
 *   exp(-\gamma\sqrt{3}\sum_{k=1}^{K}\frac{\sqrt{(x - x')^2}}{l_k}) \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @param gamma
 * @throw std::domain error if sigma <= 0, l <= 0, or x or gamma are nan of inf
 */
template <typename T_x, typename T_l, typename T_s, typename T_g>
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_l, T_s, T_g>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    gp_matern_3_2_cov(const std::vector<T_x> &x,
                      const std::vector<T_l> &length_scale,
                      const T_s &sigma = 1.0, const T_g &gamma = 1.0) {
  using std::pow;
  using std::abs;
  using std::exp;
  
  size_t x_size = x.size();
  size_t l_size = length_scale.size();

  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);

  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_not_nan("gp_matern_3_2_cov", "gamma", gamma);
  
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x", x[n]);

  Eigen::Matrix<typename stan::return_type<T_x, T_l, T_s, T_g>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3 = pow(3.0, 0.5);
  T_g neg_gamma = -1.0 * gamma;
  T_l temp;

  for (size_t i = 0; i < x_size; ++i) {
    for (size_t j = i; j < x_size; ++j) {
      temp = 0;
      for (size_t k = 0; k < l_size; ++k) {
        temp += squared_distance(x[i], x[j]) / length_scale[k];
      }
      cov(i, j) = sigma_sq * (1.0 + root_3 * gamma * pow(temp, 0.5)) *
        exp(neg_gamma * root_3 * pow(temp, 0.5));
      cov(j, i) = cov(i, j);
    }
  }
  return cov;
}

/** Returns a Matern 3/2 covariance matrix with two input vectors
 *
 * \f[ k(x, x') = \sigma^2(1 + \gamma \sqrt{3} 
 *  \frac{\sqrt{(x - x')^2}}{l}exp(-\gamma\sqrt{3}\frac{\sqrt{(x - x')^2}}{l})  \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 *
 * @param x1 std::vector of elements that can be used in stan::math::distance
 * @param x2 std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @param gamma
 * @throw std::domain error if sigma <= 0, l <= 0, or x or gamma are nan of inf
 *
 */
template <typename T_x1, typename T_x2, typename T_l, typename T_s,
          typename T_g>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_l, T_s, T_g>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_matern_3_2_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                  const T_l &length_scale, const T_s &sigma = 1.0,
                  const T_g &gamma = 1.0) {
  using std::pow;
  using std::abs;
  using std::exp;
  using stan::math::squared_distance;
  
  size_t x1_size = x1.size();
  size_t x2_size = x2.size();

  // add check same size
  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);
  
  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_not_nan("gp_matern_3_2_cov", "gamma", gamma);

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x2", x2[n]);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_g, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3_inv_l = pow(3, 0.5) / length_scale;
  T_g neg_gamma = -1.0 * gamma;

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = sigma_sq * (1.0 + root_3_inv_l * gamma * squared_distance(x1[i], x2[j])) *
        exp(neg_gamma * root_3_inv_l * squared_distance(x1[i], x2[j]));
    }
  }
  return cov;
}

/** Returns a Matern 3/2 Kernel with two input vectors with automatic 
 * relevance determination (ARD) priors
 *
 * \f[ k(x, x') = \sigma^2(1 + \gamma \sqrt{3}
 *   \frac{\sum_{k=1}^{K}\sqrt{(x - x')^2}}{l_k}
 *   exp(-\gamma\sqrt{3}\sum_{k=1}^{K}\frac{\sqrt{(x - x')^2}}{l_k}) \f]
 *
 * See Rausmussen & Williams et al 2006 Chapter 4.
 *
 * @param x1 std::vector of elements that can be used in stan::math::distance
 * @param x2 std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @param gamma
 * @throw std::domain error if sigma <= 0, l <= 0, or x or gamma are nan of inf
 *
 */
template <typename T_x1, typename T_x2, typename T_l, typename T_s,
          typename T_g>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_l, T_s, T_g>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_matern_3_2_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                  const std::vector<T_l> &length_scale, const T_s &sigma = 1.0,
                  const T_g &gamma = 1.0) {
  using std::pow;
  using std::abs;
  using std::exp;
  using stan::math::squared_distance;
  
  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  size_t l_size = length_scale.size();
  // add check same size

  check_positive("gp_matern_3_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_3_2_cov", "length-scale", length_scale);
  
  check_positive("gp_matern_3_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_3_2_cov", "marginal variance", sigma);

  check_not_nan("gp_matern_3_2_cov", "gamma", gamma);
  
  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern_3_2_cov", "x2", x2[n]);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_g, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_3 = pow(3.0, 0.5);
  T_g neg_gamma = -1.0 * gamma;
  T_l temp;

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      temp = 0;
      for (size_t k = 0; k < l_size; ++k) {
        temp += squared_distance(x1[i], x2[j]) / length_scale[k];
      }
      cov(i, j) = sigma_sq * (1.0 + root_3 * gamma * pow(temp, 0.5)) *
        exp(neg_gamma * root_3 * pow(temp, 0.5));
    }
  }
  return cov;
}
} // namespace math
} // namespace stan
#endif
