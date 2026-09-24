#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_RNG_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/binomial_rng.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a draw from a Multinomial distribution given a
 * vector of unnormalized log probabilities, a total count,
 * and a pseudo-random number generator.
 *
 * In the case of finite and -inf values this is equivalent to
 * <code>multinomial_rng(softmax(beta), N, rng)</code>.
 * Assumes a uniform distribution among all +inf values.
 *
 * @tparam RNG Type of pseudo-random number generator.
 * @param beta Vector of unnormalized log probabilities.
 * @param N Total count.
 * @param rng Pseudo-random number generator.
 * @return Multinomial random variate.
 * @throw std::domain_error if beta contains NaN
 * @throw std::domain_error if every entry of beta is negative infinity
 * @throw std::domain_error is N is less than 0.
 */
template <class RNG, typename T_beta,
          require_eigen_col_vector_t<T_beta>* = nullptr>
inline std::vector<int> multinomial_logit_rng(const T_beta& beta, int N,
                                              RNG& rng) {
  static constexpr const char* function = "multinomial_logit_rng";
  const auto& beta_ref = to_ref(beta);
  check_nonnegative(function, "number of trials variables", N);
  check_not_nan(function, "Log-probabilities parameter", beta_ref);
  if (beta_ref.size() > 0) {
    check_greater(function, "Log-probabilities parameter", beta_ref.maxCoeff(),
                  NEGATIVE_INFTY);
  }

  auto is_pos_inf = (beta_ref.array() == INFTY);
  int num_infty = is_pos_inf.count();

  plain_type_t<T_beta> theta;

  // INFTY case: uniform over the +inf entries, zero probability elsewhere
  if (num_infty > 0) {
    theta = is_pos_inf.template cast<double>() / num_infty;
  } else {
    theta = softmax(beta_ref);
  }

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
