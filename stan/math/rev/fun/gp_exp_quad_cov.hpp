#ifndef STAN_MATH_REV_FUN_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_REV_FUN_GP_EXP_QUAD_COV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <cmath>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a squared exponential kernel.
 *
 * @tparam T_x type of elements in the vector
 * @param x std::vector input that can be used in square distance
 *    Assumes each element of x is the same size
 * @param sigma standard deviation
 * @param length_scale length scale
 * @return squared distance
 * @throw std::domain_error if sigma <= 0, l <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l,
	  require_stan_scalar_t<T_sigma>* = nullptr,
	  require_stan_scalar_t<T_l>* = nullptr,
	  require_any_st_var<T_x, T_sigma, T_l>* = nullptr>
inline Eigen::Matrix<var, -1, -1>
gp_exp_quad_cov(const std::vector<T_x>& x,
		T_sigma sigma,
		T_l length_scale) {
  check_positive("gp_exp_quad_cov", "marginal standard deviation", sigma);
  check_positive("gp_exp_quad_cov", "length scale", length_scale);
  for (size_t i = 0; i < x.size(); ++i) {
    check_not_nan("gp_exp_quad_cov", "x", x[i]);
  }

  if (x.size() == 0) {
    return Eigen::Matrix<var, -1, -1>();
  }

  double l_d = value_of(length_scale);
  double sigma_d = value_of(sigma);
  double sigma_sq_d = sigma_d * sigma_d;

  arena_matrix<Eigen::MatrixXd> dist(x.size(), x.size());
  arena_matrix<Eigen::MatrixXd> res_val(x.size(), x.size());

  auto arena_x = to_arena_if<!is_constant<T_x>::value>(x);

  double inv_half_sq_l_d = 0.5 / (value_of(length_scale) * value_of(length_scale));
  for (size_t j = 0; j < x.size(); ++j) {
    for (size_t i = 0; i < j; ++i) {
      double dist_sq = squared_distance(value_of(x[i]), value_of(x[j]));
      dist(i, j) = dist_sq;
      res_val(i, j) = sigma_sq_d * std::exp(-dist_sq * inv_half_sq_l_d);
      res_val(j, i) = res_val(i, j);
      dist(j, i) = dist(i, j);
    }
  }
  for (size_t i = 0; i < x.size(); ++i) {
    dist(i, i) = 0.0;
    res_val(i, i) = sigma_sq_d;
  }

  arena_matrix<Eigen::Matrix<var,
			     Eigen::Dynamic,
			     Eigen::Dynamic>> res = res_val;

  reverse_pass_callback([=]() mutable {
    Eigen::ArrayXXd adj_times_val = res.adj().array() * res.val().array();

    /*std::cout << res.adj() << std::endl << "----" << std::endl;
    std::cout << res.val() << std::endl << "----" << std::endl;
    std::cout << adj_times_val << std::endl << "----" << std::endl;
    std::cout << dist.array() << std::endl << "----" << std::endl;*/
    
    if(!is_constant<T_x>::value)
      for(size_t i = 0; i < arena_x.size(); ++i) {
	for(size_t j = 0; j < arena_x.size(); ++j) {
	  auto adj = eval(-(value_of(arena_x[i]) - value_of(arena_x[j])) *
			  adj_times_val(i, j) / (l_d * l_d));
	  //std::cout << "(" << i << ", " << j << ") : " << adj << std::endl;
	  accumulate_adjoints(arena_x[i], adj);
	  accumulate_adjoints(arena_x[j], -adj);
	}
      }
    
    if(!is_constant<T_sigma>::value) {
      accumulate_adjoints(sigma, 2.0 * adj_times_val.sum() / sigma_d);
      /*for(size_t i = 0; i < arena_x.size(); ++i) {
	sigma
	}*/
    }

    if(!is_constant<T_l>::value)
      accumulate_adjoints(length_scale, (dist.array() * adj_times_val).sum() / (l_d * l_d * l_d));
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
