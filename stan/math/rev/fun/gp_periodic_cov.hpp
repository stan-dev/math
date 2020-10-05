#ifndef STAN_MATH_REV_FUN_GP_PERIODIC_COV_HPP
#define STAN_MATH_REV_FUN_GP_PERIODIC_COV_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/sin.hpp>
#include <stan/math/prim/fun/cos.hpp>
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
inline auto gp_periodic_cov(const std::vector<T_x>& x, T_sigma sigma, T_l l,
                            T_p p) {
  using ret_type = promote_scalar_t<var, Eigen::MatrixXd>;

  const char* fun = "gp_periodic_cov";

  check_positive(fun, "signal standard deviation", sigma);
  check_positive(fun, "length-scale", l);
  check_positive(fun, "period", p);

  for (size_t i = 0; i < x.size(); ++i) {
    check_not_nan(fun, "element of x", x[i]);
  }

  if (x.size() == 0) {
    return Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>();
  }
  double l_d = value_of(l);
  double p_d = value_of(p);
  double sigma_d = value_of(sigma);
  double sigma_squared = sigma_d * sigma_d;
  double pi_over_p = pi() / p_d;
  double negative_two_over_l_squared = -2.0 / (l_d * l_d);

  size_t N = x.size();
  size_t P = (N * N - N) / 2 + N;

  arena_t<Eigen::VectorXd> arena_dist(P);
  arena_t<Eigen::VectorXd> arena_sin_squared(P);
  arena_t<Eigen::VectorXd> arena_sin_squared_derivative(P);
  arena_t<promote_scalar_t<var, Eigen::VectorXd>> arena_res(P);

  auto arena_x = to_arena_if<!is_constant<T_x>::value>(x);

  size_t pos = 0;
  double inv_half_sq_l_d = 0.5 / (value_of(l) * value_of(l));
  for (size_t j = 0; j < N; ++j) {
    for (size_t i = 0; i <= j; ++i) {
      if (i != j) {
        if (is_stan_scalar<T_x>::value)
          arena_dist.coeffRef(pos)
              = std::abs(forward_as<double>(value_of(x[i]))
                         - forward_as<double>(value_of(x[j])));
        else if (is_eigen_col_vector<T_x>::value)
          arena_dist.coeffRef(pos)
              = forward_as<Eigen::VectorXd>(value_of(x[i]) - value_of(x[j]))
                    .norm();
        else if (is_eigen_row_vector<T_x>::value)
          arena_dist.coeffRef(pos)
              = forward_as<Eigen::RowVectorXd>(value_of(x[i]) - value_of(x[j]))
                    .norm();

        double sine = sin(pi_over_p * arena_dist.coeff(pos));
        double cosine = cos(pi_over_p * arena_dist.coeff(pos));
        double sine_squared = sine * sine;

        arena_sin_squared.coeffRef(pos) = sine_squared;

        arena_sin_squared_derivative.coeffRef(pos) = 2.0 * sine * cosine;

        arena_res.coeffRef(pos)
            = sigma_squared
              * std::exp(sine_squared * negative_two_over_l_squared);
      } else {
        arena_dist.coeffRef(pos) = 0.0;
        arena_sin_squared.coeffRef(pos) = 0.0;
        arena_sin_squared_derivative.coeffRef(pos) = 0.0;
        arena_res.coeffRef(pos) = sigma_squared;
      }

      pos++;
    }
  }

  arena_t<ret_type> res(N, N);

  pos = 0;
  for (size_t j = 0; j < N; ++j)
    for (size_t i = 0; i <= j; ++i) {
      res.coeffRef(j, i) = res.coeffRef(i, j) = arena_res.coeff(pos);
      pos++;
    }

  reverse_pass_callback([arena_res, arena_x, sigma, l, p, N, pi_over_p,
                         arena_dist, arena_sin_squared,
                         arena_sin_squared_derivative]() mutable {
    double sigma_d = value_of(sigma);
    double l_d = value_of(l);
    double p_d = value_of(p);

    double sigma_adj = 0.0;
    double l_adj = 0.0;
    double p_adj = 0.0;

    size_t pos = 0;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j <= i; ++j) {
        double adj_times_val
            = arena_res.coeff(pos).val() * arena_res.coeff(pos).adj();

        if (!is_constant<T_sigma>::value)
          sigma_adj += adj_times_val;

        if (!is_constant<T_l>::value)
          l_adj += arena_sin_squared.coeff(pos) * adj_times_val;

        if (!is_constant<T_p>::value)
          p_adj += arena_dist.coeff(pos)
                   * arena_sin_squared_derivative.coeff(pos) * adj_times_val;

        if (!is_constant<T_x>::value && i != j
            && arena_dist.coeff(pos) != 0.0) {
          auto adj = eval(-2 * pi_over_p
                          * (value_of(arena_x[i]) - value_of(arena_x[j]))
                          * arena_sin_squared_derivative(pos) * adj_times_val
                          / (arena_dist.coeff(pos) * l_d * l_d));
          forward_as<promote_scalar_t<var, T_x>>(arena_x[i]).adj() += adj;
          forward_as<promote_scalar_t<var, T_x>>(arena_x[j]).adj() -= adj;
        }
        pos++;
      }
    }

    if (!is_constant<T_sigma>::value)
      forward_as<var>(sigma).adj() += 2.0 * sigma_adj / sigma_d;

    if (!is_constant<T_l>::value)
      forward_as<var>(l).adj() += 4 * l_adj / (l_d * l_d * l_d);

    if (!is_constant<T_p>::value)
      forward_as<var>(p).adj() += 2 * pi() * p_adj / (p_d * p_d * l_d * l_d);
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
