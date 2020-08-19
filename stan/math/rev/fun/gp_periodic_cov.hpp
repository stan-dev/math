#ifndef STAN_MATH_REV_FUN_GP_PERIODIC_COV_HPP
#define STAN_MATH_REV_FUN_GP_PERIODIC_COV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <cmath>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a periodic covariance matrix \f$ \mathbf{K} \f$ using the input \f$
 * \mathbf{X} \f$. The elements of \f$ \mathbf{K} \f$ are defined as \f$
 * \mathbf{K}_{ij} = k(\mathbf{X}_i,\mathbf{X}_j), \f$ where \f$ \mathbf{X}_i
 * \f$ is the \f$i\f$-th row of \f$ \mathbf{X} \f$ and \n \f$
 * k(\mathbf{x},\mathbf{x}^\prime) = \sigma^2 \exp\left(-\frac{2\sin^2(\pi
 * |\mathbf{x}-\mathbf{x}^\prime|/p)}{\ell^2}\right), \f$ \n where \f$ \sigma^2
 * \f$, \f$ \ell \f$ and \f$ p \f$ are the signal variance, length-scale and
 * period.
 *
 * @tparam T_x type of elements in the std::vector
 * @param x std::vector of input elements.
 *   Assumes that all elements of x have the same size.
 * @param sigma standard deviation of the signal
 * @param l length-scale
 * @param p period
 * @return periodic covariance matrix
 * @throw std::domain_error if sigma <= 0, l <= 0, p <= 0, or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma, typename T_l, typename T_p,
          require_all_stan_scalar_t<T_sigma, T_l, T_p>* = nullptr,
          require_any_st_var<T_x, T_sigma, T_l, T_p>* = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
gp_periodic_cov(const std::vector<T_x> &x,
		T_sigma sigma, T_l l, T_p p) {
  const char *fun = "gp_periodic_cov";

  check_positive(fun, "signal standard deviation", sigma);
  check_positive(fun, "length-scale", l);
  check_positive(fun, "period", p);

  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i) {
    check_not_nan(fun, "element of x", x[i]);
  }

  if (x_size == 0) {
    return Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>();
  }

  double l_d = value_of(l);
  double p_d = value_of(p);
  double sigma_d = value_of(sigma);
  double sigma_squared = sigma_d * sigma_d;
  double pi_over_p = pi() / p_d;
  double negative_two_over_l_squared = -2.0 / (l_d * l_d);
  
  arena_matrix<Eigen::MatrixXd> arena_dist(x.size(), x.size());
  arena_matrix<Eigen::MatrixXd> arena_sin_squared(x.size(), x.size());
  arena_matrix<Eigen::MatrixXd> arena_sin_squared_derivative(x.size(), x.size());
  arena_matrix<Eigen::MatrixXd> res_val(x.size(), x.size());

  auto arena_x = to_arena_if<!is_constant<T_x>::value>(x);

  double inv_half_sq_l_d
      = 0.5 / (value_of(l) * value_of(l));
  for (size_t j = 0; j < x.size(); ++j) {
    for (size_t i = 0; i < j; ++i) {
      arena_dist.coeffRef(i, j) = arena_dist.coeffRef(j, i) =
	distance(value_of(x[i]), value_of(x[j]));
      arena_sin_squared.coeffRef(i, j) = arena_sin_squared.coeffRef(j, i) =
	square(sin(pi_over_p * arena_dist.coeff(i, j)));
      arena_sin_squared_derivative.coeffRef(i, j) =
	arena_sin_squared_derivative.coeffRef(j, i) = sin(2 * pi_over_p * arena_dist(i, j));

      res_val.coeffRef(i, j) = res_val.coeffRef(j, i) =
	sigma_squared * std::exp(arena_sin_squared.coeff(i, j) *
				 negative_two_over_l_squared);
    }
  }

  for (size_t i = 0; i < x.size(); ++i) {
    arena_sin_squared(i, i) = 0.0;
    arena_sin_squared_derivative(i, i) = 0.0;
    res_val(i, i) = sigma_squared;
  }

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> res
      = res_val;

  reverse_pass_callback([res, res_val, arena_x,
			 sigma, l, p,
			 sigma_d, l_d, p_d,
			 pi_over_p,
			 arena_dist,
			 arena_sin_squared,
			 arena_sin_squared_derivative]() mutable {
    Eigen::ArrayXXd adj_times_val = res.adj().array() * res_val.array();

    if (!is_constant<T_x>::value)
      for (size_t i = 0; i < arena_x.size(); ++i) {
        for (size_t j = 0; j < arena_x.size(); ++j) {
	  if(arena_dist.coeff(i, j) != 0.0) {
	    auto adj = eval(-2 * pi_over_p *
			    (value_of(arena_x[i]) - value_of(arena_x[j])) *
			    arena_sin_squared_derivative(i, j) *
			    adj_times_val(i, j) /
			    (arena_dist.coeff(i, j) * l_d * l_d));
	    using T_x_var = promote_scalar_t<var, T_x>;
	    forward_as<T_x_var>(arena_x[i]).adj() += adj;
	    forward_as<T_x_var>(arena_x[j]).adj() -= adj;
	  }
        }
      }

    if (!is_constant<T_sigma>::value)
      forward_as<var>(sigma).adj() += 2.0 * adj_times_val.sum() / sigma_d;

    if (!is_constant<T_l>::value)
      forward_as<var>(l).adj()
	+= 4 * (arena_sin_squared.array() * adj_times_val).sum() / (l_d * l_d * l_d);

    if (!is_constant<T_p>::value)
      forward_as<var>(p).adj()
	+= 2 * pi() * (arena_dist.array() * arena_sin_squared_derivative.array() * adj_times_val).sum() / (p_d * p_d * l_d * l_d);
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
