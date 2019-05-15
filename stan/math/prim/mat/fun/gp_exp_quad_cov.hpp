#ifndef STAN_MATH_PRIM_MAT_FUN_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_PRIM_MAT_FUN_GP_EXP_QUAD_COV_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/meta/get.hpp>
#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/mat/fun/divide_columns.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/mat/meta/is_constant_struct.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/err/check_nan.hpp>
#include <stan/math/opencl/gp_exp_quad_cov.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/divide_columns.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <cmath>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {
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
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    gp_exp_quad_cov(const std::vector<T_x> &x, const T_sigma &sigma,
                    const T_l &length_scale) {
  using std::exp;
  check_positive("gp_exp_quad_cov", "magnitude", sigma);
  check_positive("gp_exp_quad_cov", "length scale", length_scale);

  size_t x_size = x.size();
  Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  for (size_t n = 0; n < x.size(); ++n)
    check_not_nan("gp_exp_quad_cov", "x", x[n]);

  T_sigma sigma_sq = square(sigma);
  T_l neg_half_inv_l_sq = -0.5 / square(length_scale);

  for (size_t j = 0; j < x_size; ++j) {
    cov(j, j) = sigma_sq;
    for (size_t i = j + 1; i < x_size; ++i) {
      cov(i, j)
          = sigma_sq * exp(squared_distance(x[i], x[j]) * neg_half_inv_l_sq);
      cov(j, i) = cov(i, j);
    }
  }
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
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    gp_exp_quad_cov(const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1>> &x,
                    const T_sigma &sigma,
                    const std::vector<T_l> &length_scale) {
  using std::exp;

  check_positive_finite("gp_exp_quad_cov", "magnitude", sigma);
  check_positive_finite("gp_exp_quad_cov", "length scale", length_scale);

  size_t x_size = x.size();
  Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);
  if (x_size == 0)
    return cov;

  size_t l_size = length_scale.size();
  check_size_match("gp_exp_quad_cov", "x dimension", x[0].size(),
                   "number of length scales", l_size);

  T_sigma sigma_sq = square(sigma);
  std::vector<
      Eigen::Matrix<typename return_type<T_x, T_l>::type, Eigen::Dynamic, 1>>
      x_new = divide_columns(x, length_scale);

  for (size_t j = 0; j < x_size; ++j) {
    cov(j, j) = sigma_sq;
    for (size_t i = j + 1; i < x_size; ++i) {
      cov(i, j) = sigma_sq * exp(-0.5 * squared_distance(x_new[i], x_new[j]));
      cov(j, i) = cov(i, j);
    }
  }
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
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                const T_sigma &sigma, const T_l &length_scale) {
  using std::exp;

  const char *function_name = "gp_exp_quad_cov";
  check_positive(function_name, "magnitude", sigma);
  check_positive(function_name, "length scale", length_scale);

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);
  if (x1_size == 0 || x2_size == 0)
    return cov;

  for (size_t i = 0; i < x1_size; ++i)
    check_not_nan(function_name, "x1", x1[i]);
  for (size_t i = 0; i < x2_size; ++i)
    check_not_nan(function_name, "x2", x2[i]);

  T_sigma sigma_sq = square(sigma);
  T_l neg_half_inv_l_sq = -0.5 / square(length_scale);

  for (size_t i = 0; i < x1.size(); ++i) {
    for (size_t j = 0; j < x2.size(); ++j) {
      cov(i, j)
          = sigma_sq * exp(squared_distance(x1[i], x2[j]) * neg_half_inv_l_sq);
    }
  }
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
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_s, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<Eigen::Matrix<T_x1, Eigen::Dynamic, 1>> &x1,
                const std::vector<Eigen::Matrix<T_x2, Eigen::Dynamic, 1>> &x2,
                const T_s &sigma, const std::vector<T_l> &length_scale) {
  using std::exp;
  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  size_t l_size = length_scale.size();

  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_s, T_l>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1_size, x2_size);
  if (x1_size == 0 || x2_size == 0)
    return cov;

  const char *function_name = "gp_exp_quad_cov";
  for (size_t i = 0; i < x1_size; ++i)
    check_not_nan(function_name, "x1", x1[i]);
  for (size_t i = 0; i < x2_size; ++i)
    check_not_nan(function_name, "x2", x2[i]);
  check_positive_finite(function_name, "magnitude", sigma);
  check_positive_finite(function_name, "length scale", length_scale);
  check_size_match(function_name, "x dimension", x1[0].size(),
                   "number of length scales", l_size);
  check_size_match(function_name, "x dimension", x2[0].size(),
                   "number of length scales", l_size);

  T_s sigma_sq = square(sigma);

  std::vector<Eigen::Matrix<typename return_type<T_x1, T_l, T_s>::type,
                            Eigen::Dynamic, 1>>
      x1_new = divide_columns(x1, length_scale);
  std::vector<Eigen::Matrix<typename return_type<T_x2, T_l, T_s>::type,
                            Eigen::Dynamic, 1>>
      x2_new = divide_columns(x2, length_scale);

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = sigma_sq * exp(-0.5 * squared_distance(x1_new[i], x2_new[j]));
    }
  }
  return cov;
}

#ifdef STAN_OPENCL
/**
 * Returns a squared exponential kernel.
 *
 * @param x std::vector of scalars that can be used in square distance.
 * @param sigma marginal standard deviation or magnitude
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <>
inline Eigen::MatrixXd gp_exp_quad_cov(const std::vector<double> &x,
                                       const double &sigma,
                                       const double &length_scale) {
  const char *function_name = "gp_exp_quad_cov";
  check_positive(function_name, "magnitude", sigma);
  check_positive(function_name, "length scale", length_scale);

  const size_t x_size = x.size();
  Eigen::MatrixXd cov(x_size, x_size);
  if (x_size * x_size < opencl_context.tuning_opts().gp_exp_quad_cov_size) {
    // using explicit template args to call CPU function
    return gp_exp_quad_cov<double, double, double>(x, sigma, length_scale);
  }

  if (x_size == 0)
    return cov;

  matrix_cl x_cl(x, 1, x.size());
  check_nan(function_name, "x", x_cl);
  matrix_cl cov_cl = gp_exp_quad_cov(x_cl, sigma, length_scale);
  cov = from_matrix_cl(cov_cl);  // NOLINT

  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * @param x std::vector of scalars that can be used in square distance.
 *    This function assumes each element of x is the same size.
 * @param sigma marginal standard deviation or magnitude
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <>
inline Eigen::MatrixXd gp_exp_quad_cov(const std::vector<Eigen::VectorXd> &x,
                                       const double &sigma,
                                       const double &length_scale) {
  const char *function_name = "gp_exp_quad_cov";
  check_positive(function_name, "magnitude", sigma);
  check_positive(function_name, "length scale", length_scale);

  const size_t x_size = x.size();
  Eigen::MatrixXd cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  const size_t inner_x1_size = x[0].size();
  if (x_size * x_size * opencl_context.tuning_opts().gp_exp_quad_cov_coeff1
          + (x_size + x_size + 1) * inner_x1_size
                * opencl_context.tuning_opts().gp_exp_quad_cov_coeff2
      < 1) {
    // using explicit template args to call CPU function
    return gp_exp_quad_cov<Eigen::VectorXd, double, double>(x, sigma,
                                                            length_scale);
  }

  matrix_cl x_cl(x);
  check_nan(function_name, "x", x_cl);
  matrix_cl cov_cl = gp_exp_quad_cov(x_cl, sigma, length_scale);
  cov = from_matrix_cl(cov_cl);  // NOLINT

  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * @param x std::vector of Eigen vectors of scalars.
 * @param sigma marginal standard deviation or magnitude
 * @param length_scale std::vector length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <>
inline Eigen::MatrixXd gp_exp_quad_cov(
    const std::vector<Eigen::VectorXd> &x, const double &sigma,
    const std::vector<double> &length_scale) {
  const char *function_name = "gp_exp_quad_cov";
  check_positive_finite(function_name, "magnitude", sigma);
  check_positive_finite(function_name, "length scale", length_scale);

  const size_t x_size = x.size();
  Eigen::MatrixXd cov(x_size, x_size);

  if (x_size == 0)
    return cov;

  const size_t inner_x1_size = x[0].size();
  if (x_size * x_size * opencl_context.tuning_opts().gp_exp_quad_cov_coeff1
          + (x_size + x_size + 1) * inner_x1_size
                * opencl_context.tuning_opts().gp_exp_quad_cov_coeff2
      < 1) {
    // using explicit template args to call CPU function
    return gp_exp_quad_cov<double, double, double>(x, sigma, length_scale);
  }

  check_size_match(function_name, "x dimension", x[0].size(),
                   "number of length scales", length_scale.size());

  std::vector<Eigen::VectorXd> x_new = divide_columns(x, length_scale);

  matrix_cl x_cl(x_new);
  check_nan(function_name, "x", x_cl);
  matrix_cl cov_cl = gp_exp_quad_cov(x_cl, sigma, 1);
  cov = from_matrix_cl(cov_cl);  // NOLINT
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * This function is for the cross covariance matrix
 * needed to compute posterior predictive density.
 *
 * @param x1 std::vector of elements that can be used in square distance
 * @param x2 std::vector of elements that can be used in square distance
 * @param sigma standard deviation
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <>
inline typename Eigen::MatrixXd gp_exp_quad_cov(const std::vector<double> &x1,
                                                const std::vector<double> &x2,
                                                const double &sigma,
                                                const double &length_scale) {
  const char *function_name = "gp_exp_quad_cov";
  check_positive_finite(function_name, "magnitude", sigma);
  check_positive_finite(function_name, "length scale", length_scale);

  Eigen::MatrixXd cov(x1.size(), x2.size());
  if (cov.size() == 0)
    return cov;
  if (cov.size() < opencl_context.tuning_opts().gp_exp_quad_cov_size) {
    // using explicit template args to call CPU function
    return gp_exp_quad_cov<double, double, double, double>(x1, x2, sigma,
                                                           length_scale);
  }

  matrix_cl x1_cl(x1, 1, x1.size());
  check_nan(function_name, "x1", x1_cl);
  matrix_cl x2_cl(x2, 1, x2.size());
  check_nan(function_name, "x2", x2_cl);
  matrix_cl cov_cl = gp_exp_quad_cov(x1_cl, x2_cl, sigma, length_scale);
  cov = from_matrix_cl(cov_cl);  // NOLINT
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * This function is for the cross covariance matrix
 * needed to compute posterior predictive density.
 *
 * @param x1 std::vector of elements that can be used in square distance
 * @param x2 std::vector of elements that can be used in square distance
 *    This function assumes each element of x1 and x2 are the same size.
 * @param sigma standard deviation
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <>
inline typename Eigen::MatrixXd gp_exp_quad_cov(
    const std::vector<Eigen::VectorXd> &x1,
    const std::vector<Eigen::VectorXd> &x2, const double &sigma,
    const double &length_scale) {
  const char *function_name = "gp_exp_quad_cov";
  check_positive_finite(function_name, "magnitude", sigma);
  check_positive_finite(function_name, "length scale", length_scale);

  const size_t x1_size = x1.size();
  const size_t x2_size = x2.size();
  Eigen::MatrixXd cov(x1_size, x2_size);
  if (cov.size() == 0)
    return cov;

  const size_t inner_x1_size = x1[0].size();
  if (x1_size * x2_size * opencl_context.tuning_opts().gp_exp_quad_cov_coeff1
          + (x1_size + x2_size + 1) * inner_x1_size
                * opencl_context.tuning_opts().gp_exp_quad_cov_coeff2
      < 1) {
    // using explicit template args to call CPU function
    return gp_exp_quad_cov<Eigen::VectorXd, Eigen::VectorXd, double, double>(
        x1, x2, sigma, length_scale);
  }

  matrix_cl x1_cl(x1);
  check_nan(function_name, "x1", x1_cl);
  matrix_cl x2_cl(x2);
  check_nan(function_name, "x2", x2_cl);
  matrix_cl cov_cl = gp_exp_quad_cov(x1_cl, x2_cl, sigma, length_scale);
  cov = from_matrix_cl(cov_cl);  // NOLINT
  return cov;
}

/**
 * Returns a squared exponential kernel.
 *
 * This function is for the cross covariance
 * matrix needed to compute the posterior predictive density.
 *
 * @param x1 std::vector of Eigen vectors of scalars.
 * @param x2 std::vector of Eigen vectors of scalars.
 * @param sigma standard deviation
 * @param length_scale std::vector of length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <>
inline typename Eigen::MatrixXd gp_exp_quad_cov(
    const std::vector<Eigen::VectorXd> &x1,
    const std::vector<Eigen::VectorXd> &x2, const double &sigma,
    const std::vector<double> &length_scale) {
  const char *function_name = "gp_exp_quad_cov";
  const size_t x1_size = x1.size();
  const size_t x2_size = x2.size();
  const size_t l_size = length_scale.size();

  Eigen::MatrixXd cov(x1_size, x2_size);
  if (cov.size() == 0)
    return cov;

  const size_t inner_x1_size = x1[0].size();
  if (x1_size * x2_size * opencl_context.tuning_opts().gp_exp_quad_cov_coeff1
          + (x1_size + x2_size + 1) * inner_x1_size
                * opencl_context.tuning_opts().gp_exp_quad_cov_coeff2
      < 1) {
    // using explicit template args to call CPU function
    return gp_exp_quad_cov<double, double, double, double>(x1, x2, sigma,
                                                           length_scale);
  }

  check_positive_finite(function_name, "magnitude", sigma);
  check_positive_finite(function_name, "length scale", length_scale);
  check_size_match(function_name, "x1 dimension", x1[0].size(),
                   "number of length scales", l_size);
  check_size_match(function_name, "x2 dimension", x2[0].size(),
                   "number of length scales", l_size);

  std::vector<Eigen::VectorXd> x1_new = divide_columns(x1, length_scale);
  std::vector<Eigen::VectorXd> x2_new = divide_columns(x2, length_scale);

  matrix_cl x1_cl(x1_new);
  check_nan(function_name, "x1", x1_cl);
  matrix_cl x2_cl(x2_new);
  check_nan(function_name, "x2", x2_cl);
  matrix_cl cov_cl = gp_exp_quad_cov(x1_cl, x2_cl, sigma, 1);
  cov = from_matrix_cl(cov_cl);  // NOLINT
  return cov;
}
#endif

}  // namespace math
}  // namespace stan
#endif
