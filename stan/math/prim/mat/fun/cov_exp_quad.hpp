#ifndef STAN_MATH_PRIM_MAT_FUN_COV_EXP_QUAD_HPP
#define STAN_MATH_PRIM_MAT_FUN_COV_EXP_QUAD_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/gp_exp_quad_cov.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns a squared exponential kernel.
 *
 * @deprecated use <code>gp_exp_quad_cov</code>
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of length scale
 *
 * @param x std::vector of elements that can be used in square distance.
 *    This function assumes each element of x is the same size.
 * @param sigma standard deviation
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l>
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    cov_exp_quad(const std::vector<T_x>& x, const T_sigma& sigma,
                 const T_l& length_scale) {
  return gp_exp_quad_cov(x, sigma, length_scale);
}

/**
 * Returns a squared exponential kernel.
 *
 * @deprecated use <code>gp_exp_quad_cov</code>
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of std::vector of length scale
 *
 * @param x std::vector of elements that can be used in square distance.
 *    This function assumes each element of x is the same size.
 * @param sigma standard deviation
 * @param length_scale std::vector length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l>
inline
    typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
                           Eigen::Dynamic, Eigen::Dynamic>
    cov_exp_quad(const std::vector<T_x>& x, const T_sigma& sigma,
                 const std::vector<T_l>& length_scale) {
  return gp_exp_quad_cov(x, sigma, length_scale);
}

/**
 * Returns a squared exponential kernel.
 *
 * @deprecated use <code>gp_exp_quad_cov</code>
 *
 * @tparam T_x1 type of first std::vector of elements
 * @tparam T_x2 type of second std::vector of elements
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
cov_exp_quad(const std::vector<T_x1>& x1, const std::vector<T_x2>& x2,
             const T_sigma& sigma, const T_l& length_scale) {
  return gp_exp_quad_cov(x1, x2, sigma, length_scale);
}

/**
 * Returns a squared exponential kernel.
 *
 * @deprecated use <code>gp_exp_quad_cov</code>
 *
 * @tparam T_x1 type of first std::vector of elements
 * @tparam T_x2 type of second std::vector of elements
 * @tparam T_sigma type of sigma
 * @tparam T_l type of length scale
 *
 * @param x1 std::vector of elements that can be used in square distance
 * @param x2 std::vector of elements that can be used in square distance
 * @param sigma standard deviation
 * @param length_scale std::vector of length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<
    typename stan::return_type<T_x1, T_x2, T_sigma, T_l>::type, Eigen::Dynamic,
    Eigen::Dynamic>
cov_exp_quad(const std::vector<T_x1>& x1, const std::vector<T_x2>& x2,
             const T_sigma& sigma, const std::vector<T_l>& length_scale) {
  return gp_exp_quad_cov(x1, x2, sigma, length_scale);
}
}  // namespace math
}  // namespace stan
#endif
