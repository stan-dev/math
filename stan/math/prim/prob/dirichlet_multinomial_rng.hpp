#ifndef STAN_MATH_PRIM_PROB_DIRICHLET_MULTINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_DIRICHLET_MULTINOMIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/dirichlet_rng.hpp>
#include <stan/math/prim/prob/multinomial_rng.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a draw from a Dirichlet-Multinomial distribution with specified
 * parameters \f$\alpha\f$ and \f$N\f$ and pseudo-random number generator rng.
 *
 * The Dirichlet-Multinomial distribution is a continuous mixture of
 * Multinomial distirbutions, where the mixing distribution is the Dirichlet
 * distribution. This fact is used for generating DirMult random draws.
 * First, we sample a probability vector
 * \f$\theta \sim \mbox{Dirichlet}(\alpha)\f$.
 * Then, we sample a \f$n \sim \mbox{Multinomial}(\theta, N)\f$,
 * and return \f$n\f$.
 *
 * @tparam RNG type of pseudo-random number generator
 * @param alpha Prior sample sizes (or intensity vector).
 * @param N Number of trials.
 * @param rng Pseudo-random number generator.
 * @return A non-negative integer vector n with sum(n) = N.
 * @throw std::domain_error if any element of alpha is less than
 * or equal to 0, or infinite.
 * @throw std::domain_error is N is less than 0.
 */

template <class RNG>
inline std::vector<int> dirichlet_multinomial_rng(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, int N, RNG& rng) {
  static const char* function = "dirichlet_multinomial_rng";
  const auto& alpha_ref = to_ref(alpha);
  check_positive_finite(function, "prior size variable", alpha_ref);
  check_nonnegative(function, "number of trials variables", N);

  // special case N = 0
  if (N == 0) {
    return std::vector<int>(alpha.size(), 0);
  }

  // sample a simplex theta from the Dirichlet distribution
  auto theta = dirichlet_rng(alpha_ref, rng);

  // using the simplex theta, sample from the multinomial distribution
  return multinomial_rng(theta, N, rng);
}

}  // namespace math
}  // namespace stan
#endif
