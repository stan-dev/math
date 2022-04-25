#ifndef STAN_MATH_PRIM_PROB_INV_WISHART_CHOLESKY_RNG_HPP
#define STAN_MATH_PRIM_PROB_INV_WISHART_CHOLESKY_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/prob/wishart_cholesky_rng.hpp>
#include <stan/math/prim/prob/wishart_rng.hpp>

namespace stan {
namespace math {

template <class RNG>
inline Eigen::MatrixXd inv_wishart_cholesky_rng(double nu, const Eigen::MatrixXd& L,
                                       RNG& rng) {
  using Eigen::MatrixXd;
  static const char* function = "inv_wishart_cholesky_rng";
  index_type_t<MatrixXd> k = L.rows();
  check_greater(function, "degrees of freedom > dims - 1", nu, k - 1);
  check_cholesky_factor(function, "scale parameter", L);
  
  MatrixXd L_inv = mdivide_left_tri<Eigen::Lower>(L);
  return mdivide_left_tri<Eigen::Lower>(wishart_cholesky_rng(nu, L_inv, rng));  
}

}  // namespace math
}  // namespace stan
#endif
