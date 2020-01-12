#ifndef STAN_MATH_PRIM_PROB_INV_WISHART_RNG_HPP
#define STAN_MATH_PRIM_PROB_INV_WISHART_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inverse_spd.hpp>
#include <stan/math/prim/prob/wishart_rng.hpp>

namespace stan {
namespace math {

template <class RNG>
inline Eigen::MatrixXd inv_wishart_rng(double nu, const Eigen::MatrixXd& S,
                                       RNG& rng) {
  static const char* function = "inv_wishart_rng";

  using Eigen::MatrixXd;
  typename index_type<MatrixXd>::type k = S.rows();

  check_greater(function, "degrees of freedom > dims - 1", nu, k - 1);
  check_square(function, "scale parameter", S);

  MatrixXd S_inv = MatrixXd::Identity(k, k);
  S_inv = S.ldlt().solve(S_inv);
  MatrixXd asym = inverse_spd(wishart_rng(nu, S_inv, rng));
  return 0.5 * (asym.transpose() + asym);  // ensure symmetry
}

}  // namespace math
}  // namespace stan
#endif
