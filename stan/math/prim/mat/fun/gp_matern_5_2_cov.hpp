#ifndef STAN_MATH_PRIM_MAT_FUN_GP_MATERN_5_2_COV_HPP
#define STAN_MATH_PRIM_MAT_FUN_GP_MATERN_5_2_COV_HPP

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

template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern_5_2_cov(const std::vector<T_x> &x, const T_s &sigma,
                  const T_l &length_scale) {
  using stan::math::square;
  using stan::math::squared_distance;
  using std::exp;
  using std::pow;

  size_t x_size = x.size();
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern_5_2_cov", "x", x[n]);

  check_positive("gp_matern_5_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_5_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_5_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_5_2_cov", "length-scale", length_scale);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_5_inv_l = pow(5.0, 0.5) / length_scale;
  T_l inv_l_sq_5_3 = 5.0 / (3.0 * square(length_scale));
  T_l neg_root_5_inv_l = -1.0 * pow(5.0, 0.5) / length_scale;
  
  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_sq *
        (1.0 + root_5_inv_l * squared_distance(x[i], x[j]) +
         inv_l_sq_5_3 * pow(squared_distance(x[i], x[j]), 2.0)) *
        exp(neg_root_5_inv_l * squared_distance(x[i], x[j]));
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}

template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern_5_2_cov(const std::vector<T_x> &x, const T_s &sigma,
                  const std::vector<T_l> &length_scale) {
  using stan::math::square;
  using stan::math::squared_distance;
  using std::exp;
  using std::pow;

  size_t x_size = x.size();
  size_t l_size = length_scale.size();
  for (size_t n = 0; n < x_size; ++n)
    check_not_nan("gp_matern_5_2_cov", "x", x[n]);

  check_positive("gp_matern_5_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_5_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_5_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_5_2_cov", "length-scale", length_scale);

  Eigen::Matrix<typename stan::return_type<T_x, T_s, T_l>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_5 = pow(5.0, 0.5);
  T_l five_thirds = 5.0 / 3.0;
  T_l neg_root_5 = -1.0 * pow(5.0, 0.5);
  T_l temp;
  
  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_sq;
    for (size_t j = i + 1; j < x_size; ++j) {
      temp = 0;
      for (size_t k = 0; k < l_size; ++k)
        temp += squared_distance(x[i], x[j]) / square(length_scale[k]);
      cov(i, j) = sigma_sq *
        (1.0 + root_5 * temp + five_thirds * square(temp)) *
        exp(neg_root_5 * temp);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq;
  return cov;
}

template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x1, T_x2,
                                                         T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern_5_2_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2, 
                  const T_s &sigma, const T_l &length_scale) {
  using stan::math::square;
  using stan::math::squared_distance;
  using std::exp;
  using std::pow;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern_5_2_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern_5_2_cov", "x1", x2[n]);

  check_positive("gp_matern_5_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_5_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_5_2_cov", "length-scale", length_scale);
  check_not_nan("gp_matern_5_2_cov", "length-scale", length_scale);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_5_inv_l = pow(5.0, 0.5) / length_scale;
  T_l inv_l_sq_5_3 = 5.0 / (3.0 * square(length_scale));
  T_l neg_root_5_inv_l = -1.0 * pow(5.0, 0.5) / length_scale;
  
  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = sigma_sq *
        (1.0 + root_5_inv_l * squared_distance(x1[i], x2[j]) +
         inv_l_sq_5_3 * pow(squared_distance(x1[i], x2[j]), 2.0)) *
        exp(neg_root_5_inv_l * squared_distance(x1[i], x2[j]));
    }
  }
  return cov;
}

template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<typename stan::return_type<T_x1, T_x2,
                                                         T_s, T_l>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern_5_2_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2, 
                  const T_s &sigma, const std::vector<T_l> &length_scale) {
  using stan::math::square;
  using stan::math::squared_distance;
  using std::exp;
  using std::pow;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  size_t l_size = length_scale.size();

  for (size_t n = 0; n < x1_size; ++n)
    check_not_nan("gp_matern_5_2_cov", "x1", x1[n]);
  for (size_t n = 0; n < x2_size; ++n)
    check_not_nan("gp_matern_5_2_cov", "x1", x2[n]);

  check_positive("gp_matern_5_2_cov", "marginal variance", sigma);
  check_not_nan("gp_matern_5_2_cov", "marginal variance", sigma);

  check_positive("gp_matern_5_2_cov", "length-scale", length_scale);
  for (size_t n = 0; n < l_size; ++n)
    check_not_nan("gp_matern_5_2_cov", "length-scale", length_scale[n]);

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0)
    return cov;

  T_s sigma_sq = square(sigma);
  T_l root_5 = pow(5.0, 0.5);
  T_l five_thirds = 5.0 / 3.0;
  T_l neg_root_5 = -1.0 * pow(5.0, 0.5);
  T_l temp;
  
  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      temp = 0;
      for (size_t k = 0; k < l_size; ++k)
        temp += squared_distance(x1[i], x2[j]) / square(length_scale[k]);
      cov(i, j) = sigma_sq *
        (1.0 + root_5 * temp + five_thirds * square(temp)) *
        exp(neg_root_5 * temp);
    }
  }
  return cov;
}
}  // namespace math
}  // namespace stan
#endif
