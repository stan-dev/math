#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_RNG_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/binomial_rng.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a draw from a Multinomial distribution given a
 * a vector of unnormalized log probabilities and a pseudo-random
 * number generator.
 *
 * @tparam RNG Type of pseudo-random number generator.
 * @param beta Vector of unnormalized log probabilities.
 * @param N Total count
 * @param rng Pseudo-random number generator.
 * @return Multinomial random variate
 */
template <class RNG, typename T_beta,
          require_eigen_col_vector_t<T_beta>* = nullptr>
inline std::vector<int> multinomial_logit_rng(const T_beta& beta, int N,
                                              RNG& rng) {
  static constexpr const char* function = "multinomial_logit_rng";
  const auto& beta_ref = to_ref(beta);
  check_finite(function, "Log-probabilities parameter", beta_ref);
  check_positive(function, "number of trials variables", N);

  plain_type_t<T_beta> theta = softmax(beta_ref);
  std::vector<int> result(theta.size(), 0);
  double mass_left = 1.0;
  int n_left = N;

  for (int k = 0; n_left > 0 && k < theta.size(); ++k) {
    double p = theta.coeff(k) / mass_left;
    if (p > 1.0) {
      p = 1.0;
    }
    result[k] = binomial_rng(n_left, p, rng);
    n_left -= result[k];
    mass_left -= theta.coeff(k);
  }

  return result;
}  // namespace math

}  // namespace math
}  // namespace stan
#endif
