#ifndef STAN_MATH_PRIM_PROB_BETA_BINOMIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_BINOMIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/prob/beta_binomial_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the CDF of the Beta-Binomial distribution with given population
 * size, prior success, and prior failure parameters. Given containers of
 * matching sizes, returns the product of probabilities.
 *
 * The CDF is the exponential of `beta_binomial_lcdf`, which computes a small
 * probability from the mirrored series instead of the complement of the
 * upper tail.
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
 * @return probability or product of probabilities
 * @throw std::domain_error if N, alpha, or beta fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_N, typename T_size1, typename T_size2>
inline return_type_t<T_size1, T_size2> beta_binomial_cdf(const T_n& n,
                                                         const T_N& N,
                                                         const T_size1& alpha,
                                                         const T_size2& beta) {
  using T_N_ref = ref_type_t<T_N>;
  using T_alpha_ref = ref_type_t<T_size1>;
  using T_beta_ref = ref_type_t<T_size2>;
  static constexpr const char* function = "beta_binomial_cdf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "First prior sample size parameter", alpha,
                         "Second prior sample size parameter", beta);
  if (size_zero(n, N, alpha, beta)) {
    return 1.0;
  }

  T_N_ref N_ref = N;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_nonnegative(function, "Population size parameter", N_ref);
  check_positive_finite(function, "First prior sample size parameter",
                        alpha_ref);
  check_positive_finite(function, "Second prior sample size parameter",
                        beta_ref);

  return exp(beta_binomial_lcdf(n, N_ref, alpha_ref, beta_ref));
}

}  // namespace math
}  // namespace stan
#endif
