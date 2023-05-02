#ifndef STAN_MATH_PRIM_PROB_BINOMIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BINOMIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
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
 * Returns the log CCDF for the binomial distribution evaluated at the
 * specified success, population size, and chance of success. If given
 * containers of matching lengths, returns the log sum of probabilities.
 *
 * @tparam T_n type of successes parameter
 * @tparam T_N type of population size parameter
 * @tparam theta type of chance of success parameter
 * @param n successes parameter
 * @param N population size parameter
 * @param theta chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if N is negative
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_N, typename T_prob>
return_type_t<T_prob> binomial_lccdf(const T_n& n, const T_N& N,
                                     const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_N, T_prob>;
  using T_n_ref = ref_type_t<T_n>;
  using T_N_ref = ref_type_t<T_N>;
  using T_theta_ref = ref_type_t<T_prob>;
  using std::exp;
  using std::log;
  using std::pow;
  static const char* function = "binomial_lccdf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "Probability parameter", theta);

  T_n_ref n_ref = n;
  T_N_ref N_ref = N;
  T_theta_ref theta_ref = theta;

  check_nonnegative(function, "Population size parameter", N_ref);
  check_bounded(function, "Probability parameter", value_of(theta_ref), 0.0,
                1.0);

  if (size_zero(n, N, theta)) {
    return 0;
  }

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(theta_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_N_ref> N_vec(N_ref);
  scalar_seq_view<T_theta_ref> theta_vec(theta_ref);
  size_t max_size_seq_view = max_size(n, N, theta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined,
  // but treated as negative infinity
  for (size_t i = 0; i < stan::math::size(n); i++) {
    if (n_vec.val(i) < 0) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const T_partials_return n_dbl = n_vec.val(i);
    const T_partials_return N_dbl = N_vec.val(i);

    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (n_dbl >= N_dbl) {
      return ops_partials.build(NEGATIVE_INFTY);
    }

    const T_partials_return theta_dbl = theta_vec.val(i);
    const T_partials_return Pi
        = 1.0 - inc_beta(N_dbl - n_dbl, n_dbl + 1, 1 - theta_dbl);

    P += log(Pi);

    if (!is_constant_all<T_prob>::value) {
      const T_partials_return denom = beta(N_dbl - n_dbl, n_dbl + 1) * Pi;
      partials<0>(ops_partials)[i] += pow(theta_dbl, n_dbl)
                                      * pow(1 - theta_dbl, N_dbl - n_dbl - 1)
                                      / denom;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
