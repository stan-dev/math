#ifndef STAN_MATH_REV_FUN_GP_PERIODIC_COV_HPP
#define STAN_MATH_REV_FUN_GP_PERIODIC_COV_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <stan/math/rev/fun/square.hpp>
#include <stan/math/rev/fun/squared_distance.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/gp_exp_quad_cov.hpp>
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
template <typename T_x, typename T_sigma, require_st_arithmetic<T_x>* = nullptr,
          require_stan_scalar_t<T_sigma>* = nullptr>
inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> gp_periodic_cov(
    const std::vector<T_x>& x, const T_sigma sigma, const var l, const var p) {
  const char* fun = "gp_periodic_cov";
  check_positive(fun, "signal standard deviation", sigma);
  check_positive(fun, "length-scale", l);
  check_positive(fun, "period", p);
  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i) {
    check_not_nan(fun, "element of x", x[i]);
  }

  Eigen::Matrix<var, -1, -1> cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }
  size_t l_tri_size = x_size * (x_size - 1) / 2;
  arena_matrix<Eigen::VectorXd> dists_lin(l_tri_size);
  arena_matrix<Eigen::VectorXd> sin_dists_sq_lin(l_tri_size);
  arena_matrix<Eigen::Matrix<var, -1, 1>> cov_l_tri_lin(l_tri_size);
  arena_matrix<Eigen::Matrix<var, -1, 1>> cov_diag(
      is_constant<T_sigma>::value ? 0 : x_size);

  double sigma_sq = square(value_of(sigma));
  double pi_div_p = pi() / value_of(p);
  double neg_two_inv_l_sq = -2.0 / square(value_of(l));

  size_t block_size = 10;
  size_t pos = 0;
  for (size_t jb = 0; jb < x_size; jb += block_size) {
    size_t j_end = std::min(x_size, jb + block_size);
    size_t j_size = j_end - jb;
    cov.diagonal().segment(jb, j_size)
        = Eigen::VectorXd::Constant(j_size, sigma_sq);

    if (!is_constant<T_sigma>::value) {
      cov_diag.segment(jb, j_size) = cov.diagonal().segment(jb, j_size);
    }
    for (size_t ib = jb; ib < x_size; ib += block_size) {
      size_t i_end = std::min(x_size, ib + block_size);
      for (size_t j = jb; j < j_end; ++j) {
        for (size_t i = std::max(ib, j + 1); i < i_end; ++i) {
          double dist = distance(x[i], x[j]);
          double sin_dist = sin(pi_div_p * dist);
          double sin_dist_sq = square(sin_dist);
          dists_lin.coeffRef(pos) = dist;
          sin_dists_sq_lin.coeffRef(pos) = sin_dist_sq;
          cov_l_tri_lin.coeffRef(pos) = cov.coeffRef(j, i) = cov.coeffRef(i, j)
              = sigma_sq * exp(sin_dist_sq * neg_two_inv_l_sq);
          ++pos;
        }
      }
    }
  }

  reverse_pass_callback([cov_l_tri_lin, cov_diag, dists_lin, sin_dists_sq_lin,
                         sigma, l, p, x_size]() {
    size_t l_tri_size = x_size * (x_size - 1) / 2;
    double l_val = value_of(l);
    double p_val = value_of(p);
    double two_pi_div_p = TWO_PI / p_val;

    double adjl = 0;
    double adjsigma = 0;
    double adjp = 0;

    for (size_t pos = 0; pos < l_tri_size; ++pos) {
      double prod_add
          = cov_l_tri_lin.coeff(pos).val() * cov_l_tri_lin.coeff(pos).adj();
      adjl += prod_add * sin_dists_sq_lin.coeff(pos);
      adjsigma += prod_add;
      double dist = dists_lin.coeff(pos);
      adjp += prod_add * sin(two_pi_div_p * dist) * dist;
    }
    if (!is_constant<T_sigma>::value) {
      adjsigma += (cov_diag.val().array() * cov_diag.adj().array()).sum();
      adjoint_of(sigma) += adjsigma * 2 / value_of(sigma);
    }
    double l_sq = square(l_val);
    l.adj() += adjl * 4 / (l_sq * l_val);
    p.adj() += adjp * TWO_PI / (l_sq * square(p_val));
  });

  return cov;
}  // namespace math

}  // namespace math
}  // namespace stan
#endif
