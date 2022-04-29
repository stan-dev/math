#ifndef STAN_MATH_PRIM_PROB_INV_WISHART_CHOLESKY_RNG_HPP
#define STAN_MATH_PRIM_PROB_INV_WISHART_CHOLESKY_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/prob/wishart_cholesky_rng.hpp>
#include <stan/math/prim/prob/wishart_rng.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a random Cholesky factor of a covariance matrix
 * of the specified dimensionality drawn
 * from the inverse Wishart distribution with the specified degrees of freedom
 * using the specified random number generator.
 *
 * @tparam RNG Random number generator type
 * @param[in] nu scalar degrees of freedom
 * @param[in] L_S lower Cholesky factor of the scale matrix
 * @param[in, out] rng random-number generator
 * @return random lower Cholesky factor drawn from the given inverse Wishart
 * distribution
 * @throw std::domain_error if the scale matrix is not a Cholesky factor
 * @throw std::domain_error if the degrees of freedom is greater than k - 1
 * where k is the dimension of L
 */
template <class RNG>
inline Eigen::MatrixXd inv_wishart_cholesky_rng(double nu,
                                                const Eigen::MatrixXd& L_S,
                                                RNG& rng) {
  using Eigen::MatrixXd;
  static const char* function = "inv_wishart_cholesky_rng";
  index_type_t<MatrixXd> k = L_S.rows();
  check_greater(function, "degrees of freedom > dims - 1", nu, k - 1);
  check_positive(function, "Cholesky Scale matrix", L_S.diagonal());
  check_positive(function, "columns of Cholesky Scale matrix", L_S.cols());

  MatrixXd L_Sinv = mdivide_left_tri<Eigen::Lower>(L_S);
  return mdivide_left_tri<Eigen::Lower>(wishart_cholesky_rng(nu, L_Sinv, rng));
}

}  // namespace math
}  // namespace stan
#endif
