#ifndef STAN_MATH_PRIM_FUN_COV_EXP_QUAD_HPP
#define STAN_MATH_PRIM_FUN_COV_EXP_QUAD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/gp_exp_quad_cov.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * @deprecated use <code>gp_exp_quad_cov</code>
 */
template <typename T_x, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
cov_exp_quad(const std::vector<T_x>& x, const T_sigma& sigma,
             const T_l& length_scale) {
  return gp_exp_quad_cov(x, sigma, length_scale);
}

/**
 * @deprecated use <code>gp_exp_quad_cov</code>
 */
template <typename T_x, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x, T_sigma, T_l>, Eigen::Dynamic,
                              Eigen::Dynamic>
cov_exp_quad(const std::vector<T_x>& x, const T_sigma& sigma,
             const std::vector<T_l>& length_scale) {
  return gp_exp_quad_cov(x, sigma, length_scale);
}

/**
 * @deprecated use <code>gp_exp_quad_cov</code>
 */
template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma, T_l>,
                              Eigen::Dynamic, Eigen::Dynamic>
cov_exp_quad(const std::vector<T_x1>& x1, const std::vector<T_x2>& x2,
             const T_sigma& sigma, const T_l& length_scale) {
  return gp_exp_quad_cov(x1, x2, sigma, length_scale);
}

/**
 * @deprecated use <code>gp_exp_quad_cov</code>
 */
template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
inline typename Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma, T_l>,
                              Eigen::Dynamic, Eigen::Dynamic>
cov_exp_quad(const std::vector<T_x1>& x1, const std::vector<T_x2>& x2,
             const T_sigma& sigma, const std::vector<T_l>& length_scale) {
  return gp_exp_quad_cov(x1, x2, sigma, length_scale);
}
}  // namespace math
}  // namespace stan
#endif
