#ifndef STAN_MATH_PRIM_PROB_BINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF for the binomial distribution evaluated at the
 * specified success, population size, and chance of success. If given
 * containers of matching lengths, returns the log sum of probabilities.
 *
 * @tparam T_n type of successes parameter
 * @tparam T_N type of population size parameter
 * @tparam T_prob type of chance of success parameter
 * @param n successes parameter
 * @param N population size parameter
 * @param theta chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if n is negative or greater than N
 * @throw std::domain_error if N is negative
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_N, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_N, T_prob>* = nullptr>
return_type_t<T_prob> binomial_lpmf(const T_n& n, const T_N& N,
                                    const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_N, T_prob>;
  using T_n_ref = ref_type_t<T_n>;
  using T_N_ref = ref_type_t<T_N>;
  using T_theta_ref = ref_type_t<T_prob>;
  static const char* function = "binomial_lpmf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N,
                         "Probability parameter", theta);

  T_n_ref n_ref = n;
  T_N_ref N_ref = N;
  T_theta_ref theta_ref = theta;

  check_bounded(function, "Successes variable", value_of(n_ref), 0, N_ref);
  check_nonnegative(function, "Population size parameter", N_ref);
  check_bounded(function, "Probability parameter", value_of(theta_ref), 0.0,
                1.0);

  if (size_zero(n, N, theta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  T_partials_return logp = 0;
  auto ops_partials = make_partials_propagator(theta_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_N_ref> N_vec(N_ref);
  scalar_seq_view<T_theta_ref> theta_vec(theta_ref);
  size_t size_theta = stan::math::size(theta);
  size_t max_size_seq_view = max_size(n, N, theta);

  VectorBuilder<true, T_partials_return, T_prob> log1m_theta(size_theta);
  for (size_t i = 0; i < size_theta; ++i) {
    log1m_theta[i] = log1m(theta_vec.val(i));
  }

  if (include_summand<propto>::value) {
    for (size_t i = 0; i < max_size_seq_view; ++i) {
      logp += binomial_coefficient_log(N_vec[i], n_vec[i]);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; ++i) {
    if (N_vec[i] != 0) {
      if (n_vec[i] == 0) {
        logp += N_vec[i] * log1m_theta[i];
      } else if (n_vec[i] == N_vec[i]) {
        logp += n_vec[i] * log(theta_vec.val(i));
      } else {
        logp += n_vec[i] * log(theta_vec.val(i))
                + (N_vec[i] - n_vec[i]) * log1m_theta[i];
      }
    }
  }

  if (!is_constant_all<T_prob>::value) {
    if (size_theta == 1) {
      T_partials_return sum_n = 0;
      T_partials_return sum_N = 0;
      for (size_t i = 0; i < max_size_seq_view; ++i) {
        sum_n += n_vec[i];
        sum_N += N_vec[i];
      }
      const T_partials_return theta_dbl = theta_vec.val(0);
      if (sum_N != 0) {
        if (sum_n == 0) {
          partials<0>(ops_partials)[0] -= sum_N / (1.0 - theta_dbl);
        } else if (sum_n == sum_N) {
          partials<0>(ops_partials)[0] += sum_n / theta_dbl;
        } else {
          partials<0>(ops_partials)[0]
              += sum_n / theta_dbl - (sum_N - sum_n) / (1.0 - theta_dbl);
        }
      }
    } else {
      for (size_t i = 0; i < max_size_seq_view; ++i) {
        const T_partials_return theta_dbl = theta_vec.val(i);
        if (N_vec[i] != 0) {
          if (n_vec[i] == 0) {
            partials<0>(ops_partials)[i] -= N_vec[i] / (1.0 - theta_dbl);
          } else if (n_vec[i] == N_vec[i]) {
            partials<0>(ops_partials)[i] += n_vec[i] / theta_dbl;
          } else {
            partials<0>(ops_partials)[i]
                += n_vec[i] / theta_dbl
                   - (N_vec[i] - n_vec[i]) / (1.0 - theta_dbl);
          }
        }
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
