#ifndef STAN_MATH_PRIM_PROB_LKJ_CORR_RNG_HPP
#define STAN_MATH_PRIM_PROB_LKJ_CORR_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/prim/prob/lkj_corr_cholesky_rng.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a random correlation matrix (symmetric, positive
 * definite, unit diagonal) of the specified dimensionality drawn
 * from the LKJ distribution with the specified degrees of freedom
 * using the specified random number generator.
 *
 * @tparam RNG Random number generator type.
 * @param[in] K Number of rows and columns of generated matrix.
 * @param[in] eta Degrees of freedom for LKJ distribution.
 * @param[in, out] rng Random-number generator to use.
 * @return Random variate with specified distribution.
 * @throw std::domain_error If the shape parameter is not positive.
 */
template <class RNG>
inline Eigen::MatrixXd lkj_corr_rng(size_t K, double eta, RNG& rng) {
  static constexpr const char* function = "lkj_corr_rng";
  check_positive(function, "Shape parameter", eta);
  return multiply_lower_tri_self_transpose(lkj_corr_cholesky_rng(K, eta, rng));
}

}  // namespace math
}  // namespace stan
#endif
