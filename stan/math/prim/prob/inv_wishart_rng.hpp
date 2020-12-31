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
  using Eigen::MatrixXd;
  static const char* function = "inv_wishart_rng";
  index_type_t<MatrixXd> k = S.rows();
  check_greater(function, "degrees of freedom > dims - 1", nu, k - 1);
  check_square(function, "scale parameter", S);
  check_symmetric(function, "scale parameter", S);

  MatrixXd S_inv = MatrixXd::Identity(k, k);
  Eigen::LDLT<Eigen::MatrixXd> ldlt_of_S = S.ldlt();
  check_pos_definite("inv_wishart_rng", "scale parameter", ldlt_of_S);
  S_inv = ldlt_of_S.solve(S_inv);
  MatrixXd asym = inverse_spd(wishart_rng(nu, S_inv, rng));
  return 0.5 * (asym.transpose() + asym);  // ensure symmetry
}

}  // namespace math
}  // namespace stan
#endif
