#ifndef STAN_MATH_REV_FUN_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_REV_FUN_GP_EXP_QUAD_COV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
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
template <typename T_x, typename T_sigma, require_st_arithmetic<T_x>* = nullptr,
          require_stan_scalar_t<T_sigma>* = nullptr>
inline Eigen::Matrix<var, -1, -1> gp_exp_quad_cov(const std::vector<T_x>& x,
                                                  const T_sigma sigma,
                                                  const var length_scale) {
  check_positive("gp_exp_quad_cov", "sigma", sigma);
  check_positive("gp_exp_quad_cov", "length_scale", length_scale);
  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i) {
    check_not_nan("gp_exp_quad_cov", "x", x[i]);
  }

  Eigen::Matrix<var, -1, -1> cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }
  size_t l_tri_size = x_size * (x_size - 1) / 2;
  arena_matrix<Eigen::VectorXd> sq_dists_lin(l_tri_size);
  arena_matrix<Eigen::Matrix<var, -1, 1>> cov_l_tri_lin(l_tri_size);
  arena_matrix<Eigen::Matrix<var, -1, 1>> cov_diag(
      is_constant<T_sigma>::value ? 0 : x_size);

  double l_val = value_of(length_scale);
  double sigma_sq = square(value_of(sigma));
  double neg_half_inv_l_sq = -0.5 / square(l_val);

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
          sq_dists_lin.coeffRef(pos) = squared_distance(x[i], x[j]);
          cov_l_tri_lin.coeffRef(pos) = cov.coeffRef(j, i) = cov.coeffRef(i, j)
              = sigma_sq * exp(sq_dists_lin.coeff(pos) * neg_half_inv_l_sq);
          pos++;
        }
      }
    }
  }

  reverse_pass_callback(
      [cov_l_tri_lin, cov_diag, sq_dists_lin, sigma, length_scale, x_size]() {
        size_t l_tri_size = x_size * (x_size - 1) / 2;
        double adjl = 0;
        double adjsigma = 0;
        for (Eigen::Index pos = 0; pos < l_tri_size; pos++) {
          double prod_add
              = cov_l_tri_lin.coeff(pos).val() * cov_l_tri_lin.coeff(pos).adj();
          adjl += prod_add * sq_dists_lin.coeff(pos);
          if (!is_constant<T_sigma>::value) {
            adjsigma += prod_add;
          }
        }
        if (!is_constant<T_sigma>::value) {
          adjsigma += (cov_diag.val().array() * cov_diag.adj().array()).sum();
          adjoint_of(sigma) += adjsigma * 2 / value_of(sigma);
        }
        double l_val = value_of(length_scale);
        length_scale.adj() += adjl / (l_val * l_val * l_val);
      });

  return cov;
}

}  // namespace math
}  // namespace stan
#endif
