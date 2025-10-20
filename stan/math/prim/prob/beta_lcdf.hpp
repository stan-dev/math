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
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

#include <cmath>

namespace stan {
namespace math {

/**
 * Beta log CDF.
 *
 * @tparam T_y type of y
 * @tparam T_scale_succ type of success parameter
 * @tparam T_scale_fail type of failure parameter
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
inline return_type_t<T_y, T_scale_succ, T_scale_fail> beta_lcdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta_param) {
  using T_partials_return = partials_return_t<T_y, T_scale_succ, T_scale_fail>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale_succ>;
  using T_beta_ref = ref_type_t<T_scale_fail>;
  static constexpr const char* function = "beta_lcdf";

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

  const size_t size_alpha = stan::math::size(alpha);
  const size_t size_beta = stan::math::size(beta_param);
  const size_t size_alpha_beta = max_size(alpha, beta_param);
  const size_t N = max_size(y, alpha, beta_param);

  // Allocate digamma buffers only if alpha/beta contain any autodiff scalars.
  constexpr bool need_digamma = is_any_autodiff_v<T_scale_succ, T_scale_fail>;
  VectorBuilder<need_digamma, T_partials_return, T_scale_succ> digamma_alpha(
      size_alpha);
  VectorBuilder<need_digamma, T_partials_return, T_scale_fail> digamma_beta(
      size_beta);
  VectorBuilder<need_digamma, T_partials_return, T_scale_succ, T_scale_fail>
      digamma_sum(size_alpha_beta);

  if constexpr (is_any_autodiff_v<T_scale_succ, T_scale_fail>) {
    for (size_t i = 0; i < size_alpha; ++i) {
      digamma_alpha[i] = digamma(alpha_vec.val(i));
    }
    for (size_t i = 0; i < size_beta; ++i) {
      digamma_beta[i] = digamma(beta_vec.val(i));
    }
    for (size_t i = 0; i < size_alpha_beta; ++i) {
      digamma_sum[i] = digamma(alpha_vec.val(i) + beta_vec.val(i));
    }
  }

  for (size_t n = 0; n < N; ++n) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return alpha_dbl = alpha_vec.val(n);
    const T_partials_return beta_dbl = beta_vec.val(n);
    const T_partials_return betafunc_dbl = beta(alpha_dbl, beta_dbl);
    const T_partials_return Pn = inc_beta(alpha_dbl, beta_dbl, y_dbl);

    const T_partials_return inv_Pn
        = is_any_autodiff_v<T_y, T_scale_succ, T_scale_fail> ? inv(Pn) : 0;

    cdf_log += log(Pn);

    if constexpr (is_any_autodiff_v<T_y>) {
      partials<0>(ops_partials)[n] += pow(1 - y_dbl, beta_dbl - 1)
                                      * pow(y_dbl, alpha_dbl - 1) * inv_Pn
                                      / betafunc_dbl;
    }

    if constexpr (is_any_autodiff_v<T_scale_succ, T_scale_fail>) {
      T_partials_return g1 = 0;
      T_partials_return g2 = 0;
      grad_reg_inc_beta(g1, g2, alpha_dbl, beta_dbl, y_dbl, digamma_alpha[n],
                        digamma_beta[n], digamma_sum[n], betafunc_dbl);
      if constexpr (is_any_autodiff_v<T_scale_succ>) {
        partials<1>(ops_partials)[n] += g1 * inv_Pn;
      }
      if constexpr (is_any_autodiff_v<T_scale_fail>) {
        partials<2>(ops_partials)[n] += g2 * inv_Pn;
      }
    }
  }

  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif  // STAN_MATH_PRIM_PROB_BETA_LCDF_HPP
