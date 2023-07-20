#ifndef STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/hypergeometric_3F2.hpp>
#include <stan/math/prim/fun/grad_F32.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CCDF of the Beta-Binomial distribution with given population
 * size, prior success, and prior failure parameters. Given containers of
 * matching sizes, returns the log sum of probabilities.
 *
 * @tparam T_n type of success parameter
 * @tparam T_N type of population size parameter
 * @tparam T_size1 type of prior success parameter
 * @tparam T_size2 type of prior failure parameter
 *
 * @param n success parameter
 * @param N population size parameter
 * @param alpha prior success parameter
 * @param beta prior failure parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if N, alpha, or beta fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_N, typename T_size1, typename T_size2>
return_type_t<T_size1, T_size2> beta_binomial_lccdf(const T_n& n, const T_N& N,
                                                    const T_size1& alpha,
                                                    const T_size2& beta) {
  using T_partials_return = partials_return_t<T_n, T_N, T_size1, T_size2>;
  using std::exp;
  using std::log;
  using T_N_ref = ref_type_t<T_N>;
  using T_alpha_ref = ref_type_t<T_size1>;
  using T_beta_ref = ref_type_t<T_size2>;
  static constexpr const char* function = "beta_binomial_lccdf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "First prior sample size parameter", alpha,
                         "Second prior sample size parameter", beta);
  if (size_zero(n, N, alpha, beta)) {
    return 0;
  }

  T_N_ref N_ref = N;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_nonnegative(function, "Population size parameter", N_ref);
  check_positive_finite(function, "First prior sample size parameter",
                        alpha_ref);
  check_positive_finite(function, "Second prior sample size parameter",
                        beta_ref);

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(alpha_ref, beta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_N_ref> N_vec(N_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t max_size_seq_view = max_size(n, N, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(n); i++) {
    if (n_vec.val(i) < 0) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as neg infinity
    if (n_vec.val(i) >= N_vec.val(i)) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return n_dbl = n_vec.val(i);
    const T_partials_return N_dbl = N_vec.val(i);
    const T_partials_return alpha_dbl = alpha_vec.val(i);
    const T_partials_return beta_dbl = beta_vec.val(i);
    const T_partials_return mu = alpha_dbl + n_dbl + 1;
    const T_partials_return nu = beta_dbl + N_dbl - n_dbl - 1;
    const T_partials_return one = 1;

    const T_partials_return F = hypergeometric_3F2(
        {one, mu, -N_dbl + n_dbl + 1}, {n_dbl + 2, 1 - nu}, one);
    T_partials_return C = lbeta(nu, mu) - lbeta(alpha_dbl, beta_dbl)
                          - lbeta(N_dbl - n_dbl, n_dbl + 2);
    C = F * exp(C) / (N_dbl + 1);

    const T_partials_return Pi = C;

    P += log(Pi);

    T_partials_return digammaDiff
        = is_constant_all<T_size1, T_size2>::value
              ? 0
              : digamma(alpha_dbl + beta_dbl) - digamma(mu + nu);

    T_partials_return dF[6];
    if (!is_constant_all<T_size1, T_size2>::value) {
      grad_F32(dF, one, mu, -N_dbl + n_dbl + 1, n_dbl + 2, 1 - nu, one);
    }
    if (!is_constant_all<T_size1>::value) {
      partials<0>(ops_partials)[i]
          += digamma(mu) - digamma(alpha_dbl) + digammaDiff + dF[1] / F;
    }
    if (!is_constant_all<T_size2>::value) {
      partials<1>(ops_partials)[i]
          += digamma(nu) - digamma(beta_dbl) + digammaDiff - dF[4] / F;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
