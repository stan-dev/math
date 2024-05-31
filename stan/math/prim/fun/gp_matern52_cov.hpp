#ifndef STAN_MATH_PRIM_FUN_GP_MATERN52_COV_HPP
#define STAN_MATH_PRIM_FUN_GP_MATERN52_COV_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/divide_columns.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a Matern 5/2 covariance matrix with one input vector
 *
 * \f[ k(x, x') = \sigma^2\bigg(1 +
 *  \frac{\sqrt{5}d(x, x')}{l} + \frac{5d(x, x')^2}{3l^2}\bigg)
 *   exp\bigg(-\frac{5 d(x, x')}{l}\bigg)
 * \f]
 *
 * where \f$ d(x, x') \f$ is the Euclidean distance.
 *
 * @tparam T_x type of elements contained in vector x
 * @tparam T_s type of element of sigma, the magnitude
 * @tparam T_l type of elements of length scale
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
gp_matern52_cov(const std::vector<T_x> &x, const T_s &sigma,
                const T_l &length_scale) {
  using std::exp;
  using std::sqrt;

  size_t x_size = x.size();
  Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0) {
    return cov;
  }

  for (size_t n = 0; n < x_size; ++n) {
    check_not_nan("gp_matern52_cov", "x", x[n]);
  }

  check_positive_finite("gp_matern52_cov", "magnitude", sigma);
  check_positive_finite("gp_matern52_cov", "length scale", length_scale);

  T_s sigma_sq = square(sigma);
  T_l root_5_inv_l = sqrt(5.0) / length_scale;
  T_l inv_l_sq_5_3 = 5.0 / (3.0 * square(length_scale));

  size_t block_size = 10;
  for (size_t jb = 0; jb < x_size; jb += block_size) {
    for (size_t ib = jb; ib < x_size; ib += block_size) {
      size_t j_end = std::min(x_size, jb + block_size);
      for (size_t j = jb; j < j_end; ++j) {
        cov.coeffRef(j, j) = sigma_sq;
        size_t i_end = std::min(x_size, ib + block_size);
        for (size_t i = std::max(ib, j + 1); i < i_end; ++i) {
          return_type_t<T_x> sq_distance = squared_distance(x[i], x[j]);
          return_type_t<T_x> dist = sqrt(sq_distance);
          cov.coeffRef(j, i) = cov.coeffRef(i, j)
              = sigma_sq
                * (1.0 + root_5_inv_l * dist + inv_l_sq_5_3 * sq_distance)
                * exp(-root_5_inv_l * dist);
        }
      }
    }
  }
  return cov;
}

/**
 * Returns a Matern 5/2 covariance matrix with one input vector
 * with automatic relevance determination (ARD).
 *
 * \f[ k(x, x') = \sigma^2\bigg(1 +
 *  \sqrt{5}\sqrt{\sum_{k=1}^{K}\frac{d(x, x')^2}{l_k^2}} +
 *  \frac{5}{3} \sqrt{\sum_{k=1}^{K}\frac{d(x, x')^2}{l_k^2}}\bigg)
 *  exp\bigg(-\frac{5}{3}\bigg(\sqrt{\sum_{k=1}^K{\frac{d(x, x')^2}{l_k^2}}
 *  }\bigg)\bigg)
 * \f]
 *
 * where \f$ d(x, x') \f$ is the Euclidean distance.
 *
 * @tparam T_x type of elements contained in vector x
 * @tparam T_s type of element of sigma, the magnitude
 * @tparam T_l type of elements in vector of length scale
 *
 * @param x std::vector of elements that can be used in stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 * @throw std::invalid_argument if length scale size != dimension of x
 */
template <typename T_x, typename T_s, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
gp_matern52_cov(const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1>> &x,
                const T_s &sigma, const std::vector<T_l> &length_scale) {
  using std::exp;

  size_t x_size = x.size();
  Eigen::Matrix<return_type_t<T_x, T_s, T_l>, Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }

  size_t l_size = length_scale.size();
  for (size_t n = 0; n < x_size; ++n) {
    check_not_nan("gp_matern52_cov", "x", x[n]);
  }

  check_positive_finite("gp_matern52_cov", "magnitude", sigma);
  check_positive_finite("gp_matern52_cov", "length scale", length_scale);
  check_size_match("gp_matern52_cov", "x dimension", x[0].size(),
                   "number of length scales", l_size);

  T_s sigma_sq = square(sigma);
  double root_5 = sqrt(5.0);
  double five_thirds = 5.0 / 3.0;

  std::vector<Eigen::Matrix<return_type_t<T_x, T_l>, -1, 1>> x_new
      = divide_columns(x, length_scale);
  size_t block_size = 10;

  for (size_t jb = 0; jb < x_size; jb += block_size) {
    for (size_t ib = jb; ib < x_size; ib += block_size) {
      size_t j_end = std::min(x_size, jb + block_size);
      for (size_t j = jb; j < j_end; ++j) {
        cov.coeffRef(j, j) = sigma_sq;
        size_t i_end = std::min(x_size, ib + block_size);
        for (size_t i = std::max(ib, j + 1); i < i_end; ++i) {
          return_type_t<T_x, T_l> sq_distance
              = squared_distance(x_new[i], x_new[j]);
          return_type_t<T_x, T_l> dist = sqrt(sq_distance);
          cov.coeffRef(j, i) = cov.coeffRef(i, j)
              = sigma_sq * (1.0 + root_5 * dist + five_thirds * sq_distance)
                * exp(-root_5 * dist);
        }
      }
    }
  }
  return cov;
}

/**
 * Returns a Matern 5/2 covariance matrix with two different input vectors
 *
 * \f[ k(x, x') = \sigma^2\bigg(1 +
 *  \frac{\sqrt{5}d(x, x')}{l} + \frac{5d(x, x')^2}{3l^2}\bigg)
 *   exp\bigg(-\frac{5 d(x, x')}{l}\bigg)
 * \f]
 *
 * where \f$ d(x, x') \f$ is the Euclidean distance.
 *
 * @tparam T_x1 type of elements contained in vector x1
 * @tparam T_x2 type of elements contained in vector x2
 * @tparam T_s type of element of sigma, the magnitude
 * @tparam T_l type of elements of length scale
 *
 * @param x1 std::vector of elements that can be used in
 * stan::math::distance
 * @param x2 std::vector of elements that can be used in
 * stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 *
 */
template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x1, T_x2, T_s, T_l>,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern52_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                const T_s &sigma, const T_l &length_scale) {
  using std::exp;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  Eigen::Matrix<return_type_t<T_x1, T_x2, T_s, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  for (size_t n = 0; n < x1_size; ++n) {
    check_not_nan("gp_matern52_cov", "x1", x1[n]);
  }
  for (size_t n = 0; n < x2_size; ++n) {
    check_not_nan("gp_matern52_cov", "x1", x2[n]);
  }

  check_positive_finite("gp_matern52_cov", "magnitude", sigma);
  check_positive_finite("gp_matern52_cov", "length scale", length_scale);

  T_s sigma_sq = square(sigma);
  T_l root_5_inv_l = sqrt(5.0) / length_scale;
  T_l inv_l_sq_5_3 = 5.0 / (3.0 * square(length_scale));

  size_t block_size = 10;

  for (size_t ib = 0; ib < x1_size; ib += block_size) {
    for (size_t jb = 0; jb < x2_size; jb += block_size) {
      size_t j_end = std::min(x2_size, jb + block_size);
      for (size_t j = jb; j < j_end; ++j) {
        size_t i_end = std::min(x1_size, ib + block_size);
        for (size_t i = ib; i < i_end; ++i) {
          return_type_t<T_x1, T_x2> sq_distance
              = squared_distance(x1[i], x2[j]);
          return_type_t<T_x1, T_x2> dist = sqrt(sq_distance);
          cov.coeffRef(i, j)
              = sigma_sq
                * (1.0 + root_5_inv_l * dist + inv_l_sq_5_3 * sq_distance)
                * exp(-root_5_inv_l * dist);
        }
      }
    }
  }
  return cov;
}

/**
 * Returns a Matern 5/2 covariance matrix with two input vectors
 * with automatic relevance determination (ARD).
 *
 * \f[ k(x, x') = \sigma^2\bigg(1 +
 *  \sqrt{5}\sqrt{\sum_{k=1}^{K}\frac{d(x, x')^2}{l_k^2}} +
 *  \frac{5}{3} \sqrt{\sum_{k=1}^{K}\frac{d(x, x')^2}{l_k^2}}\bigg)
 *  exp\bigg(-\frac{5}{3}\bigg(\sqrt{\sum_{k=1}^K{\frac{d(x, x')^2}{l_k^2}}
 *  }\bigg)\bigg)
 * \f]
 *
 * where \f$ d(x, x') \f$ is the Euclidean distance.
 *
 * @tparam T_x1 type of elements contained in vector x1
 * @tparam T_x2 type of elements contained in vector x2
 * @tparam T_s type of element of sigma, the  magnitude
 * @tparam T_l type of elements in vector of length scale
 *
 * @param x1 std::vector of elements that can be used in
 * stan::math::distance
 * @param x2 std::vector of elements that can be used in
 * stan::math::distance
 * @param length_scale length scale
 * @param sigma standard deviation that can be used in stan::math::square
 * @throw std::domain error if sigma <= 0, l <= 0, or x is nan or inf
 * @throw std::invalid_argument if length scale size != dimension of x1 or
 * x2
 *
 */
template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x1, T_x2, T_s, T_l>,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_matern52_cov(const std::vector<Eigen::Matrix<T_x1, Eigen::Dynamic, 1>> &x1,
                const std::vector<Eigen::Matrix<T_x2, Eigen::Dynamic, 1>> &x2,
                const T_s &sigma, const std::vector<T_l> &length_scale) {
  using std::exp;

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  Eigen::Matrix<return_type_t<T_x1, T_x2, T_s, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  size_t l_size = length_scale.size();

  for (size_t n = 0; n < x1_size; ++n) {
    check_not_nan("gp_matern52_cov", "x1", x1[n]);
  }
  for (size_t n = 0; n < x2_size; ++n) {
    check_not_nan("gp_matern52_cov", "x1", x2[n]);
  }

  check_positive_finite("gp_matern52_cov", "magnitude", sigma);
  check_positive_finite("gp_matern52_cov", "length scale", length_scale);
  check_size_match("gp_matern52_cov", "x dimension", x1[0].size(),
                   "number of length scales", l_size);
  check_size_match("gp_matern52_cov", "x dimension", x2[0].size(),
                   "number of length scales", l_size);

  T_s sigma_sq = square(sigma);
  double root_5 = sqrt(5.0);
  double five_thirds = 5.0 / 3.0;

  std::vector<Eigen::Matrix<return_type_t<T_x1, T_l>, -1, 1>> x1_new
      = divide_columns(x1, length_scale);
  std::vector<Eigen::Matrix<return_type_t<T_x2, T_l>, -1, 1>> x2_new
      = divide_columns(x2, length_scale);
  size_t block_size = 10;

  for (size_t ib = 0; ib < x1_size; ib += block_size) {
    for (size_t jb = 0; jb < x2_size; jb += block_size) {
      size_t j_end = std::min(x2_size, jb + block_size);
      for (size_t j = jb; j < j_end; ++j) {
        size_t i_end = std::min(x1_size, ib + block_size);
        for (size_t i = ib; i < i_end; ++i) {
          return_type_t<T_x1, T_x2, T_l> sq_distance
              = squared_distance(x1_new[i], x2_new[j]);
          return_type_t<T_x1, T_x2, T_l> dist = sqrt(sq_distance);
          cov.coeffRef(i, j)
              = sigma_sq * (1.0 + root_5 * dist + five_thirds * sq_distance)
                * exp(-root_5 * dist);
        }
      }
    }
  }
  return cov;
}
}  // namespace math
}  // namespace stan
#endif
