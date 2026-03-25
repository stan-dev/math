#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_CDF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>
#include <stan/math/prim/prob/neg_binomial_cdf.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/subtract.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the CDF of the geometric distribution. Given containers of
 * matching sizes, returns the product of probabilities.
 *
 * Delegates to the negative binomial CDF with alpha = 1 and
 * beta = theta / (1 - theta).
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if theta is not in (0, 1]
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_prob>
inline return_type_t<T_prob> geometric_cdf(const T_n& n,
                                           const T_prob& theta) {
  using T_n_ref = ref_type_t<T_n>;
  using T_prob_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_cdf";

  check_consistent_sizes(function, "Outcome variable", n,
                         "Success probability parameter", theta);
  if (size_zero(n, theta)) {
    return 1.0;
  }

  T_n_ref n_ref = n;
  T_prob_ref theta_ref = theta;
  check_bounded(function, "Success probability parameter",
                value_of(theta_ref), 0.0, 1.0);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  for (int i = 0; i < stan::math::size(n); i++) {
    if (n_vec.val(i) < 0) {
      return 0.0;
    }
  }

  // theta = 1 => CDF is always 1 for n >= 0
  scalar_seq_view<T_prob_ref> theta_vec(theta_ref);
  bool all_theta_one = true;
  for (size_t i = 0; i < stan::math::size(theta); i++) {
    if (value_of(theta_vec[i]) != 1.0) {
      all_theta_one = false;
      break;
    }
  }
  if (all_theta_one) {
    return 1.0;
  }

  if constexpr (is_stan_scalar_v<T_prob>) {
    const auto beta = theta_ref / (1.0 - theta_ref);
    return neg_binomial_cdf(n_ref, 1, beta);
  } else if constexpr (is_std_vector_v<T_prob>) {
    std::vector<value_type_t<T_prob>> beta;
    beta.reserve(stan::math::size(theta));
    for (size_t i = 0; i < stan::math::size(theta); i++) {
      beta.push_back(theta_vec[i] / (1.0 - theta_vec[i]));
    }
    return neg_binomial_cdf(n_ref, 1, beta);
  } else {
    const auto beta = elt_divide(theta_ref, subtract(1.0, theta_ref));
    return neg_binomial_cdf(n_ref, 1, beta);
  }
}

}  // namespace math
}  // namespace stan
#endif
