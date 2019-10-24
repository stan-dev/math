#ifndef STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

/**
 * Returns the log PMF for the binomial distribution evaluated at the
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
 * @throw std::domain_error if n is negative or greater than N
 * @throw std::domain_error if N is negative
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_N, typename T_prob>
return_type_t<T_prob> binomial_lpmf(const T_n& n, const T_N& N,
                                    const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_N, T_prob>;

  static const char* function = "binomial_lpmf";

  if (size_zero(n, N, theta)) {
    return 0.0;
  }

  T_partials_return logp = 0;
  check_bounded(function, "Successes variable", n, 0, N);
  check_nonnegative(function, "Population size parameter", N);
  check_finite(function, "Probability parameter", theta);
  check_bounded(function, "Probability parameter", theta, 0.0, 1.0);
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "Probability parameter", theta);

  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_N> N_vec(N);
  scalar_seq_view<T_prob> theta_vec(theta);
  size_t size = max_size(n, N, theta);

  operands_and_partials<T_prob> ops_partials(theta);

  if (include_summand<propto>::value) {
    for (size_t i = 0; i < size; ++i) {
      logp += binomial_coefficient_log(N_vec[i], n_vec[i]);
    }
  }

  VectorBuilder<true, T_partials_return, T_prob> log1m_theta(length(theta));
  for (size_t i = 0; i < length(theta); ++i) {
    log1m_theta[i] = log1m(value_of(theta_vec[i]));
  }

  for (size_t i = 0; i < size; ++i) {
    logp += multiply_log(n_vec[i], value_of(theta_vec[i]))
            + (N_vec[i] - n_vec[i]) * log1m_theta[i];
  }

  if (length(theta) == 1) {
    T_partials_return temp1 = 0;
    T_partials_return temp2 = 0;
    for (size_t i = 0; i < size; ++i) {
      temp1 += n_vec[i];
      temp2 += N_vec[i] - n_vec[i];
    }
    if (!is_constant_all<T_prob>::value) {
      ops_partials.edge1_.partials_[0]
          += temp1 / value_of(theta_vec[0])
             - temp2 / (1.0 - value_of(theta_vec[0]));
    }
  } else {
    if (!is_constant_all<T_prob>::value) {
      for (size_t i = 0; i < size; ++i) {
        ops_partials.edge1_.partials_[i]
            += n_vec[i] / value_of(theta_vec[i])
               - (N_vec[i] - n_vec[i]) / (1.0 - value_of(theta_vec[i]));
      }
    }
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_N, typename T_prob>
inline return_type_t<T_prob> binomial_lpmf(const T_n& n, const T_N& N,
                                           const T_prob& theta) {
  return binomial_lpmf<false>(n, N, theta);
}

}  // namespace math
}  // namespace stan
#endif
