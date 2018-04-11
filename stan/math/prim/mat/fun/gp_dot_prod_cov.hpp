#ifndef STAN_MATH_PRIM_MAT_FUN_COV_DOT_PROD_HPP
#define STAN_MATH_PRIM_MAT_FUN_COV_DOT_PROD_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/fun/dot_product.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a dot product kernel.
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma type of sigma
 *
 * @param x std::vector of elements that can be used in dot product
 *    This function assumes each element of x is the same size.
 * @param sigma
 * @return dot product kernel
 * @throw std::domain_error if sigma < 0, nan, inf or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma>
inline typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<T_x> &x, const T_sigma &sigma) {
  using stan::math::dot_product;
  using stan::math::square;

  check_not_nan("gp_dot_prod_cov", "sigma", sigma);
  check_nonnegative("gp_dot_prod_cov", "sigma", sigma);
  check_finite("gp_dot_prod_cov", "sigma", sigma);

  for (size_t n = 0; n < x.size(); ++n) {
    check_not_nan("gp_dot_prod_cov", "x", x[n]);
    check_finite("gp_dot_prod_cov", "x", x[n]);
  }

  Eigen::Matrix<typename stan::return_type<T_x, T_sigma>::type, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x.size(), x.size());
  int x_size = x.size();
  if (x_size == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);

  for (int i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_sq + dot_product(x[i], x[i]);
    for (int j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_sq + dot_product(x[i], x[j]);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1)
      = sigma_sq + dot_product(x[x_size - 1], x[x_size - 1]);
  return cov;
}

/**
 * Returns a dot product kernel.
 *
 * @tparam T_x type of std::vector of double
 * @tparam T_sigma type of sigma
 *
 * @param x std::vector of elements that can be used in transpose
 *   and multiply
 *    This function assumes each element of x is the same size.
 * @param sigma
 * @return dot product kernel
 * @throw std::domain_error if sigma < 0, nan, inf or
 *   x is nan or infinite
 */
template <typename T_sigma>
inline typename Eigen::Matrix<typename stan::return_type<double, T_sigma>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<double> &x, const T_sigma &sigma) {
  using stan::math::dot_product;
  using stan::math::square;

  check_not_nan("gp_dot_prod_cov", "sigma", sigma);
  check_nonnegative("gp_dot_prod_cov", "sigma", sigma);
  check_finite("gp_dot_prod_cov", "sigma", sigma);

  for (size_t n = 0; n < x.size(); ++n) {
    check_not_nan("gp_dot_prod_cov", "x", x[n]);
    check_finite("gp_dot_prod_cov", "x", x[n]);
  }

  Eigen::Matrix<typename stan::return_type<double, T_sigma>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x.size(), x.size());
  int x_size = x.size();
  if (x_size == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);

  for (int i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_sq + x[i] * x[i];
    for (int j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_sq + x[i] * x[j];
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_sq + x[x_size - 1] * x[x_size - 1];
  return cov;
}

/**
 * Returns a dot product kernel.
 *
 * @tparam T_x1 type of first std::vector of elements
 * @tparam T_x2 type of second std::vector of elements
 * @tparam T_sigma type of sigma
 *
 * @param x1 std::vector of elements that can be used in dot_product
 * @param x2 std::vector of elements that can be used in dot_product
 * @param sigma
 * @return dot product kernel
 * @throw std::domain_error if sigma < 0, nan or inf
 *   or if x1 or x2 are nan or inf
 */
template <typename T_x1, typename T_x2, typename T_sigma>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_sigma>::type, Eigen::Dynamic,
    Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                const T_sigma &sigma) {
  using stan::math::dot_product;
  using stan::math::square;

  check_not_nan("gp_dot_prod_cov", "sigma", sigma);
  check_nonnegative("gp_dot_prod_cov", "sigma", sigma);
  check_finite("gp_dot_prod_cov", "sigma", sigma);

  for (size_t n = 0; n < x1.size(); ++n) {
    check_not_nan("gp_dot_prod_cov", "x1", x1[n]);
    check_finite("gp_dot_prod_cov", "x1", x1[n]);
  }
  for (size_t n = 0; n < x2.size(); ++n) {
    check_not_nan("gp_dot_prod_cov", "x2", x2[n]);
    check_finite("gp_dot_prod_cov", "x2", x2[n]);
  }
  Eigen::Matrix<typename stan::return_type<T_x1, T_x2, T_sigma>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1.size(), x2.size());

  if (x1.size() == 0 || x2.size() == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);

  for (size_t i = 0; i < x1.size(); ++i) {
    for (size_t j = 0; j < x2.size(); ++j) {
      cov(i, j) = sigma_sq + dot_product(x1[i], x2[j]);
    }
  }
  return cov;
}

/**
 * Returns a dot product kernel.
 *
 * @tparam T_x1 type of first std::vector of double
 * @tparam T_x2 type of second std::vector of double
 * @tparam T_sigma type of sigma
 *
 * @param x1 std::vector of elements that can be used in dot_product
 * @param x2 std::vector of elements that can be used in dot_product
 * @param sigma
 * @return dot product kernel
 * @throw std::domain_error if sigma < 0, nan or inf
 *   or if x1 or x2 are nan or inf
 */
template <typename T_sigma>
inline typename Eigen::Matrix<typename stan::return_type<double, T_sigma>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<double> &x1, const std::vector<double> &x2,
                const T_sigma &sigma) {
  using stan::math::square;

  check_not_nan("gp_dot_prod_cov", "sigma", sigma);
  check_nonnegative("gp_dot_prod_cov", "sigma", sigma);
  check_finite("gp_dot_prod_cov", "sigma", sigma);

  for (size_t n = 0; n < x1.size(); ++n) {
    check_not_nan("gp_dot_prod_cov", "x1", x1[n]);
    check_finite("gp_dot_prod_cov", "x1", x1[n]);
  }
  for (size_t n = 0; n < x2.size(); ++n) {
    check_not_nan("gp_dot_prod_cov", "x2", x2[n]);
    check_finite("gp_dot_prod_cov", "x2", x2[n]);
  }
  Eigen::Matrix<typename stan::return_type<double, T_sigma>::type,
                Eigen::Dynamic, Eigen::Dynamic>
      cov(x1.size(), x2.size());

  if (x1.size() == 0 || x2.size() == 0)
    return cov;

  T_sigma sigma_sq = square(sigma);

  for (size_t i = 0; i < x1.size(); ++i) {
    for (size_t j = 0; j < x2.size(); ++j) {
      cov(i, j) = sigma_sq + x1[i] * x2[j];
    }
  }
  return cov;
}
}  // namespace math
}  // namespace stan
#endif
