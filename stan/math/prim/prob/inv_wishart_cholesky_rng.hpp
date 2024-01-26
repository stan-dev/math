#ifndef STAN_MATH_PRIM_PROB_INV_WISHART_CHOLESKY_RNG_HPP
#define STAN_MATH_PRIM_PROB_INV_WISHART_CHOLESKY_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_right_tri_low.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a random Cholesky factor of a covariance matrix
 * of the specified dimensionality drawn
 * from the inverse Wishart distribution with the specified degrees of freedom
 * using the specified random number generator.
 *
 * Axen, Seth D. "Efficiently generating inverse-Wishart matrices and their
 * Cholesky factors." arXiv preprint arXiv:2310.15884 (2023).
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
  check_square(function, "Cholesky Scale matrix", L_S);
  check_greater(function, "degrees of freedom > dims - 1", nu, k - 1);
  check_cholesky_factor(function, "Cholesky Scale matrix", L_S);

  MatrixXd B = MatrixXd::Zero(k, k);
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < j; ++i) {
      B(j, i) = normal_rng(0, 1, rng);
    }
    B(j, j) = std::sqrt(chi_square_rng(nu - k + j + 1, rng));
  }

  return mdivide_right_tri_low(L_S, B);
}

}  // namespace math
}  // namespace stan
#endif
