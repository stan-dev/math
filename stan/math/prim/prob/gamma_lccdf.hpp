#ifndef STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fma.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/grad_reg_lower_inc_gamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/fun/log_gamma_q_dgamma.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <optional>

namespace stan {
namespace math {
namespace internal {

/**
 * Computes log q and d(log q) / d(alpha) using continued fraction.
 */
template <bool any_fvar, bool partials_fvar, typename T_shape, typename T1, typename T2>
inline std::optional<std::pair<return_type_t<T1, T2>, return_type_t<T1, T2>>>
eval_q_cf(const T1& alpha, const T2& beta_y) {
  using scalar_t = return_type_t<T1, T2>;
  using ret_t = std::pair<scalar_t, scalar_t>;
  if constexpr (!any_fvar && is_autodiff_v<T_shape>) {
    std::pair<double, double> log_q_result
        = log_gamma_q_dgamma(value_of_rec(alpha), value_of_rec(beta_y));
    if (likely(std::isfinite(value_of_rec(log_q_result.first)))) {
      return std::optional{log_q_result};
    } else {
      return std::optional<ret_t>{std::nullopt};
    }
  } else {
    ret_t out{internal::log_q_gamma_cf(alpha, beta_y), 0.0};
    if (unlikely(!std::isfinite(value_of_rec(out.first)))) {
      return std::optional<ret_t>{std::nullopt};
    }
    if constexpr (is_autodiff_v<T_shape>) {
      if constexpr (!partials_fvar) {
        out.second
            = grad_reg_inc_gamma(alpha, beta_y, tgamma(alpha), digamma(alpha))
              / exp(out.first);
      } else {
        auto alpha_unit = alpha;
        alpha_unit.d_ = 1;
        auto beta_y_unit = beta_y;
        beta_y_unit.d_ = 0;
        auto log_Q_fvar = internal::log_q_gamma_cf(alpha_unit, beta_y_unit);
        out.second = log_Q_fvar.d_;
      }
    }
    return std::optional{out};
  }
}

/**
 * Computes log q and d(log q) / d(alpha) using log1m.
 */
template <bool partials_fvar, typename T_shape, typename T1, typename T2>
inline std::optional<std::pair<return_type_t<T1, T2>, return_type_t<T1, T2>>>
eval_q_log1m(const T1& alpha, const T2& beta_y) {
  using scalar_t = return_type_t<T1, T2>;
  using ret_t = std::pair<scalar_t, scalar_t>;
  ret_t out{log1m(gamma_p(alpha, beta_y)), 0.0};
  if (unlikely(!std::isfinite(value_of_rec(out.first)))) {
    return std::optional<ret_t>{std::nullopt};
  }
  if constexpr (is_autodiff_v<T_shape>) {
    if constexpr (partials_fvar) {
      auto alpha_unit = alpha;
      alpha_unit.d_ = 1;
      auto beta_unit = beta_y;
      beta_unit.d_ = 0;
      auto log_Q_fvar = log1m(gamma_p(alpha_unit, beta_unit));
      out.second = log_Q_fvar.d_;
    } else {
      out.second
          = -grad_reg_lower_inc_gamma(alpha, beta_y) / exp(out.first);
    }
  }
  return std::optional{out};
}
}  // namespace internal

template <typename T_y, typename T_shape, typename T_inv_scale>
inline return_type_t<T_y, T_shape, T_inv_scale> gamma_lccdf(
    const T_y& y, const T_shape& alpha, const T_inv_scale& beta) {
  using std::exp;
  using std::log;
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;
  static constexpr const char* function = "gamma_lccdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);
  check_nonnegative(function, "Random variable", y_ref);

  if (size_zero(y, alpha, beta)) {
    return 0;
  }

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  const size_t N = max_size(y, alpha, beta);

  constexpr bool any_fvar = is_fvar<scalar_type_t<T_y>>::value
                            || is_fvar<scalar_type_t<T_shape>>::value
                            || is_fvar<scalar_type_t<T_inv_scale>>::value;
  constexpr bool partials_fvar = is_fvar<T_partials_return>::value;

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    const T_partials_return y_val = y_vec.val(n);
    if (y_val == 0.0) {
      continue;
    }
    if (y_val == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return alpha_val = alpha_vec.val(n);
    const T_partials_return beta_val = beta_vec.val(n);

    const T_partials_return beta_y = beta_val * y_val;
    if (beta_y == INFTY) {
      return ops_partials.build(negative_infinity());
    }
    std::optional<std::pair<T_partials_return, T_partials_return>> result;
    if (beta_y > alpha_val + 1.0) {
      result = internal::eval_q_cf<any_fvar, partials_fvar, T_shape>(alpha_val, beta_y);
    } else {
      result = internal::eval_q_log1m<partials_fvar, T_shape>(alpha_val, beta_y);
      if (!result && beta_y > 0.0) {
        // Fallback to continued fraction if log1m fails
        result = internal::eval_q_cf<any_fvar, partials_fvar, T_shape>(alpha_val, beta_y);
      }
    }
    if (unlikely(!result)) {
      return ops_partials.build(negative_infinity());
    }

    P += result->first;

    if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_inv_scale>) {
      const T_partials_return log_y = log(y_val);
      const T_partials_return alpha_minus_one = fma(alpha_val, log_y, -log_y);

      const T_partials_return log_pdf = alpha_val * log(beta_val)
                                        - lgamma(alpha_val) + alpha_minus_one
                                        - beta_y;

      const T_partials_return hazard = exp(log_pdf - result->first);  // f/Q

      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] -= hazard;
      }
      if constexpr (is_autodiff_v<T_inv_scale>) {
        partials<2>(ops_partials)[n] -= (y_val / beta_val) * hazard;
      }
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<1>(ops_partials)[n] += result->second;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan

#endif
