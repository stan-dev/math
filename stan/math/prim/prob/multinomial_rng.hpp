#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/binomial_rng.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a draw from a Multinomial distribution given a
 * probability simplex, a total count, and a pseudo-random
 * number generator.
 *
 * @tparam RNG Type of pseudo-random number generator.
 * @param theta Vector of normalized probabilities.
 * @param N Total count.
 * @param rng Pseudo-random number generator.
 * @return Multinomial random variate.
 * @throw std::domain_error if theta is not a simplex.
 * @throw std::domain_error is N is less than 0.
 */
template <class T_theta, class RNG,
          require_eigen_col_vector_t<T_theta>* = nullptr>
inline std::vector<int> multinomial_rng(const T_theta& theta, int N, RNG& rng) {
  static constexpr const char* function = "multinomial_rng";
  const auto& theta_ref = to_ref(theta);
  check_simplex(function, "Probabilities parameter", theta_ref);
  check_nonnegative(function, "number of trials variables", N);

  std::vector<int> result(theta.size(), 0);
  double mass_left = 1.0;
  int n_left = N;
  for (int k = 0; n_left > 0 && k < theta.size(); ++k) {
    double p = theta_ref.coeff(k) / mass_left;
    if (p > 1.0) {
      p = 1.0;
    }
    result[k] = binomial_rng(n_left, p, rng);
    n_left -= result[k];
    mass_left -= theta_ref(k);
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
