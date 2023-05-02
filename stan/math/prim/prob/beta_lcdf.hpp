#ifndef STAN_MATH_PRIM_PROB_BETA_LCDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the beta log cumulative distribution function for the given
 * probability, success, and failure parameters.  Any arguments other
 * than scalars must be containers of the same size.  With non-scalar
 * arguments, the return is the sum of the log cdfs with scalars
 * broadcast as necessary.
 *
 * @tparam T_y type of y
 * @tparam T_scale_succ type of success parameter
 * @tparam T_scale_fail type of failure parameter
 * @param y (Sequence of) scalar(s) between zero and one
 * @param alpha (Sequence of) success parameter(s)
 * @param beta_param (Sequence of) failure parameter(s)
 * @return log probability or sum of log of probabilities
 * @throw std::domain_error if alpha or beta is negative
 * @throw std::domain_error if y is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
return_type_t<T_y, T_scale_succ, T_scale_fail> beta_lcdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta_param) {
  using T_partials_return = partials_return_t<T_y, T_scale_succ, T_scale_fail>;
  using std::exp;
  using std::log;
  using std::pow;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale_succ>;
  using T_beta_ref = ref_type_t<T_scale_fail>;
  static const char* function = "beta_lcdf";
  check_consistent_sizes(function, "Random variable", y,
                         "First shape parameter", alpha,
                         "Second shape parameter", beta_param);
  if (size_zero(y, alpha, beta_param)) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta_param;
  check_positive_finite(function, "First shape parameter", alpha_ref);
  check_positive_finite(function, "Second shape parameter", beta_ref);
  check_bounded(function, "Random variable", value_of(y_ref), 0, 1);

  T_partials_return cdf_log(0.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t size_alpha = stan::math::size(alpha);
  size_t size_beta = stan::math::size(beta_param);
  size_t size_alpha_beta = max_size(alpha, beta_param);
  size_t N = max_size(y, alpha, beta_param);

  VectorBuilder<!is_constant_all<T_scale_succ, T_scale_fail>::value,
                T_partials_return, T_scale_succ>
      digamma_alpha(size_alpha);
  VectorBuilder<!is_constant_all<T_scale_succ, T_scale_fail>::value,
                T_partials_return, T_scale_fail>
      digamma_beta(size_beta);
  VectorBuilder<!is_constant_all<T_scale_succ, T_scale_fail>::value,
                T_partials_return, T_scale_succ, T_scale_fail>
      digamma_sum(size_alpha_beta);

  if (!is_constant_all<T_scale_succ, T_scale_fail>::value) {
    for (size_t i = 0; i < size_alpha; i++) {
      digamma_alpha[i] = digamma(alpha_vec.val(i));
    }
    for (size_t i = 0; i < size_beta; i++) {
      digamma_beta[i] = digamma(beta_vec.val(i));
    }
    for (size_t i = 0; i < size_alpha_beta; i++) {
      digamma_sum[i] = digamma(alpha_vec.val(i) + beta_vec.val(i));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return alpha_dbl = alpha_vec.val(n);
    const T_partials_return beta_dbl = beta_vec.val(n);
    const T_partials_return betafunc_dbl = beta(alpha_dbl, beta_dbl);
    const T_partials_return Pn = inc_beta(alpha_dbl, beta_dbl, y_dbl);
    const T_partials_return inv_Pn
        = is_constant_all<T_y, T_scale_succ, T_scale_fail>::value ? 0 : inv(Pn);

    cdf_log += log(Pn);

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[n] += pow(1 - y_dbl, beta_dbl - 1)
                                      * pow(y_dbl, alpha_dbl - 1) * inv_Pn
                                      / betafunc_dbl;
    }

    T_partials_return g1 = 0;
    T_partials_return g2 = 0;

    if (!is_constant_all<T_scale_succ, T_scale_fail>::value) {
      grad_reg_inc_beta(g1, g2, alpha_dbl, beta_dbl, y_dbl, digamma_alpha[n],
                        digamma_beta[n], digamma_sum[n], betafunc_dbl);
    }
    if (!is_constant_all<T_scale_succ>::value) {
      partials<1>(ops_partials)[n] += g1 * inv_Pn;
    }
    if (!is_constant_all<T_scale_fail>::value) {
      partials<2>(ops_partials)[n] += g2 * inv_Pn;
    }
  }

  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
