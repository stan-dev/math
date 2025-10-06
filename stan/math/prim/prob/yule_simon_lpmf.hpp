#ifndef STAN_MATH_PRIM_PROB_YULE_SIMON_LPMF_HPP
#define STAN_MATH_PRIM_PROB_YULE_SIMON_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the Yule-Simon distribution with shape parameter.
 * Given containers of matching sizes, returns the log sum of probabilities.
 *
 * @tparam T_n type of outcome variable
 * @tparam T_alpha type of shape parameter
 *
 * @param n outcome variable
 * @param alpha shape parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if alpha fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_alpha,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_alpha> * = nullptr>
inline return_type_t<T_alpha> yule_simon_lpmf(const T_n &n,
                                              const T_alpha &alpha) {
  using std::log;
  using T_partials_return = partials_return_t<T_n, T_alpha>;
  using T_n_ref = ref_type_t<T_n>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  static constexpr const char *function = "yule_simon_lpmf";
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha);
  if (size_zero(n, alpha)) {
    return 0.0;
  }

  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;
  check_greater_or_equal(function, "Outcome variable", n_ref, 1);
  check_positive_finite(function, "Shape parameter", alpha_ref);

  if constexpr (!include_summand<propto, T_alpha>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(alpha_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  const size_t max_size_seq_view = max_size(n_ref, alpha_ref);
  T_partials_return logp(0.0);
  if constexpr (include_summand<propto>::value) {
    if constexpr (is_stan_scalar_v<T_n>) {
      logp += lgamma(n_ref) * max_size_seq_view;
    }
  }
  for (size_t i = 0; i < max_size_seq_view; i++) {
    if constexpr (include_summand<propto>::value) {
      if constexpr (!is_stan_scalar_v<T_n>) {
        logp += lgamma(n_vec.val(i));
      }
    }
    T_partials_return alpha_plus_one = alpha_vec.val(i) + 1.0;
    logp += log(alpha_vec.val(i)) + lgamma(alpha_plus_one)
            - lgamma(n_vec.val(i) + alpha_plus_one);
    if constexpr (is_autodiff_v<T_alpha>) {
      partials<0>(ops_partials)[i] += 1.0 / alpha_vec.val(i)
                                      + digamma(alpha_plus_one)
                                      - digamma(n_vec.val(i) + alpha_plus_one);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_alpha>
inline return_type_t<T_alpha> yule_simon_lpmf(const T_n &n,
                                              const T_alpha &alpha) {
  return yule_simon_lpmf<false>(n, alpha);
}

}  // namespace math
}  // namespace stan
#endif
