#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_LPMF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/prob/neg_binomial_lpmf.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/subtract.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the geometric distribution. If containers
 * of matching sizes are supplied, returns the log sum of probabilities.
 *
 * The geometric distribution is a special case of the negative
 * binomial distribution with alpha = 1 and beta = theta / (1 - theta).
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not in (0, 1]
 * @throw std::domain_error if n is negative
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
inline return_type_t<T_prob> geometric_lpmf(const T_n& n,
                                            const T_prob& theta) {
  using T_n_ref = ref_type_t<T_n>;
  using T_prob_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_lpmf";
  check_consistent_sizes(function, "Outcome variable", n,
                         "Success probability parameter", theta);
  if (size_zero(n, theta)) {
    return 0.0;
  }

  T_n_ref n_ref = n;
  T_prob_ref theta_ref = theta;
  check_nonnegative(function, "Outcome variable", n_ref);
  check_bounded(function, "Success probability parameter",
                value_of(theta_ref), 0.0, 1.0);

  // theta = 1 => deterministic: P(0) = 1, P(n>0) = 0
  // Cannot delegate since beta = theta / (1 - theta) diverges
  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_prob_ref> theta_vec(theta_ref);
  const size_t max_sz = max_size(n_ref, theta_ref);
  for (size_t i = 0; i < max_sz; i++) {
    if (value_of(theta_vec[i]) == 1.0) {
      if (n_vec[i] > 0) {
        return negative_infinity();
      }
    }
  }
  bool all_theta_one = true;
  for (size_t i = 0; i < stan::math::size(theta); i++) {
    if (value_of(theta_vec[i]) != 1.0) {
      all_theta_one = false;
      break;
    }
  }
  if (all_theta_one) {
    return 0.0;
  }

  // geometric(theta) = neg_binomial(1, theta / (1 - theta))
  if constexpr (is_stan_scalar_v<T_prob>) {
    const auto beta = theta_ref / (1.0 - theta_ref);
    return neg_binomial_lpmf<propto>(n_ref, 1, beta);
  } else if constexpr (is_std_vector_v<T_prob>) {
    std::vector<value_type_t<T_prob>> beta;
    beta.reserve(stan::math::size(theta));
    for (size_t i = 0; i < stan::math::size(theta); i++) {
      beta.push_back(theta_vec[i] / (1.0 - theta_vec[i]));
    }
    return neg_binomial_lpmf<propto>(n_ref, 1, beta);
  } else {
    const auto beta = elt_divide(theta_ref, subtract(1.0, theta_ref));
    return neg_binomial_lpmf<propto>(n_ref, 1, beta);
  }
}

template <typename T_n, typename T_prob>
inline return_type_t<T_prob> geometric_lpmf(const T_n& n,
                                            const T_prob& theta) {
  return geometric_lpmf<false>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
