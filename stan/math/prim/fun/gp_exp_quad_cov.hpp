#ifndef STAN_MATH_PRIM_FUN_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_PRIM_FUN_GP_EXP_QUAD_COV_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/divide_columns.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

namespace internal {

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x type for each scalar
 * @tparam T_sigma type of parameter sigma
 * @tparam T_l type of parameter length scale
 *
 * @param x std::vector of scalars that can be used in square distance.
 *    This function assumes each element of x is the same size.
 * @param sigma_sq square root of the marginal standard deviation or magnitude
 * @param neg_half_inv_l_sq The half negative inverse of the length scale
 * @return squared distance
 */
template <typename T_x, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<T_x> &x, const T_sigma &sigma_sq,
                const T_l &neg_half_inv_l_sq) {
  using std::exp;
  const size_t x_size = x.size();
  Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);
  cov.diagonal().array() = sigma_sq;
  size_t block_size = 10;
  for (size_t jb = 0; jb < x.size(); jb += block_size) {
    for (size_t ib = jb; ib < x.size(); ib += block_size) {
      size_t j_end = std::min(x_size, jb + block_size);
      for (size_t j = jb; j < j_end; ++j) {
        size_t i_end = std::min(x_size, ib + block_size);
        for (size_t i = std::max(ib, j + 1); i < i_end; ++i) {
          cov.coeffRef(i, j)
              = sigma_sq
                * exp(squared_distance(x[i], x[j]) * neg_half_inv_l_sq);
        }
      }
    }
  }
  cov.template triangularView<Eigen::Upper>() = cov.transpose();
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * This function is for the cross covariance matrix
 * needed to compute posterior predictive density.
 *
 * @tparam T_x1 type of first std::vector of scalars
 * @tparam T_x2 type of second std::vector of scalars
 *    This function assumes each element of x1 and x2 are the same size.
 * @tparam T_sigma type of sigma
 * @tparam T_l type of of length scale
 *
 * @param x1 std::vector of elements that can be used in square distance
 * @param x2 std::vector of elements that can be used in square distance
 * @param sigma_sq square root of the marginal standard deviation or magnitude
 * @param neg_half_inv_l_sq The half negative inverse of the length scale
 * @return squared distance
 */
template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma, T_l>,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                const T_sigma &sigma_sq, const T_l &neg_half_inv_l_sq) {
  using std::exp;
  Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1.size(), x2.size());
  size_t block_size = 10;

  for (size_t ib = 0; ib < x1.size(); ib += block_size) {
    for (size_t jb = 0; jb < x2.size(); jb += block_size) {
      size_t j_end = std::min(x2.size(), jb + block_size);
      for (size_t j = jb; j < j_end; ++j) {
        size_t i_end = std::min(x1.size(), ib + block_size);
        for (size_t i = ib; i < i_end; ++i) {
          cov.coeffRef(i, j)
              = sigma_sq
                * exp(squared_distance(x1[i], x2[j]) * neg_half_inv_l_sq);
        }
      }
    }
  }
  return cov;
}
}  // namespace internal

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x type for each scalar
 * @tparam T_sigma type of parameter sigma
 * @tparam T_l type of parameter length scale
 *
 * @param x std::vector of scalars that can be used in square distance.
 *    This function assumes each element of x is the same size.
 * @param sigma marginal standard deviation or magnitude
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<T_x> &x, const T_sigma &sigma,
                const T_l &length_scale) {
  check_positive("gp_exp_quad_cov", "magnitude", sigma);
  check_positive("gp_exp_quad_cov", "length scale", length_scale);

  const size_t x_size = x.size();
  Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0) {
    return cov;
  }

  for (size_t n = 0; n < x_size; ++n) {
    check_not_nan("gp_exp_quad_cov", "x", x[n]);
  }

  cov = internal::gp_exp_quad_cov(x, square(sigma),
                                  -0.5 / square(length_scale));
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x type for each scalar
 * @tparam T_sigma type of parameter sigma
 * @tparam T_l type of each length scale parameter
 *
 * @param x std::vector of Eigen vectors of scalars.
 * @param sigma marginal standard deviation or magnitude
 * @param length_scale std::vector length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<Eigen::Matrix<T_x, -1, 1>> &x,
                const T_sigma &sigma, const std::vector<T_l> &length_scale) {
  check_positive_finite("gp_exp_quad_cov", "magnitude", sigma);
  check_positive_finite("gp_exp_quad_cov", "length scale", length_scale);

  size_t x_size = x.size();
  Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }

  check_size_match("gp_exp_quad_cov", "x dimension", x[0].size(),
                   "number of length scales", length_scale.size());
  cov = internal::gp_exp_quad_cov(divide_columns(x, length_scale),
                                  square(sigma), -0.5);
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * This function is for the cross covariance matrix
 * needed to compute posterior predictive density.
 *
 * @tparam T_x1 type of first std::vector of scalars
 * @tparam T_x2 type of second std::vector of scalars
 *    This function assumes each element of x1 and x2 are the same size.
 * @tparam T_sigma type of sigma
 * @tparam T_l type of of length scale
 *
 * @param x1 std::vector of elements that can be used in square distance
 * @param x2 std::vector of elements that can be used in square distance
 * @param sigma standard deviation
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma, T_l>,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                const T_sigma &sigma, const T_l &length_scale) {
  const char *function_name = "gp_exp_quad_cov";
  check_positive(function_name, "magnitude", sigma);
  check_positive(function_name, "length scale", length_scale);

  const size_t x1_size = x1.size();
  const size_t x2_size = x2.size();
  Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);
  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  for (size_t i = 0; i < x1_size; ++i) {
    check_not_nan(function_name, "x1", x1[i]);
  }
  for (size_t i = 0; i < x2_size; ++i) {
    check_not_nan(function_name, "x2", x2[i]);
  }

  cov = internal::gp_exp_quad_cov(x1, x2, square(sigma),
                                  -0.5 / square(length_scale));
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * This function is for the cross covariance
 * matrix needed to compute the posterior predictive density.
 *
 * @tparam T_x1 type of first std::vector of elements
 * @tparam T_x2 type of second std::vector of elements
 * @tparam T_s type of sigma
 * @tparam T_l type of length scale
 *
 * @param x1 std::vector of Eigen vectors of scalars.
 * @param x2 std::vector of Eigen vectors of scalars.
 * @param sigma standard deviation
 * @param length_scale std::vector of length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x1, typename T_x2, typename T_s, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x1, T_x2, T_s, T_l>,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<Eigen::Matrix<T_x1, -1, 1>> &x1,
                const std::vector<Eigen::Matrix<T_x2, -1, 1>> &x2,
                const T_s &sigma, const std::vector<T_l> &length_scale) {
  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  size_t l_size = length_scale.size();

  Eigen::Matrix<return_type_t<T_x1, T_x2, T_s, T_l>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  const char *function_name = "gp_exp_quad_cov";
  for (size_t i = 0; i < x1_size; ++i) {
    check_not_nan(function_name, "x1", x1[i]);
  }
  for (size_t i = 0; i < x2_size; ++i) {
    check_not_nan(function_name, "x2", x2[i]);
  }
  check_positive_finite(function_name, "magnitude", sigma);
  check_positive_finite(function_name, "length scale", length_scale);
  check_size_match(function_name, "x dimension", x1[0].size(),
                   "number of length scales", l_size);
  check_size_match(function_name, "x dimension", x2[0].size(),
                   "number of length scales", l_size);
  cov = internal::gp_exp_quad_cov(divide_columns(x1, length_scale),
                                  divide_columns(x2, length_scale),
                                  square(sigma), -0.5);
  return cov;
}

}  // namespace math
}  // namespace stan
#endif
