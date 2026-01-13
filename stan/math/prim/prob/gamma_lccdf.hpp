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

namespace stan {
namespace math {
namespace internal {
  template <typename T>
  struct Q_eval {
    T log_Q{0.0};
    T dlogQ_dalpha{0.0};
    bool ok{false};
  };

  /**
   * Computes log q and d(log q) / d(alpha) using continued fraction.
   */
  template <typename T, typename T_shape,
            bool any_fvar, bool partials_fvar>
  static inline Q_eval<T> eval_q_cf(const T& alpha,
                                    const T& beta_y) {
    Q_eval<T> out;
    if constexpr (!any_fvar && is_autodiff_v<T_shape>) {
      auto log_q_result = log_gamma_q_dgamma(value_of_rec(alpha), value_of_rec(beta_y));
      out.log_Q = log_q_result.log_q;
      out.dlogQ_dalpha = log_q_result.dlog_q_da;
    } else {
      out.log_Q = internal::log_q_gamma_cf(alpha, beta_y);
      if constexpr (is_autodiff_v<T_shape>) {
        if constexpr (!partials_fvar) {
          out.dlogQ_dalpha
            = grad_reg_inc_gamma(alpha, beta_y, tgamma(alpha),
                                 digamma(alpha)) / exp(out.log_Q);
        } else {
          T alpha_unit = alpha;
          alpha_unit.d_ = 1;
          T beta_y_unit = beta_y;
          beta_y_unit.d_ = 0;
          T log_Q_fvar = internal::log_q_gamma_cf(alpha_unit, beta_y_unit);
          out.dlogQ_dalpha = log_Q_fvar.d_;
        }
      }
    }

    out.ok = std::isfinite(value_of_rec(out.log_Q));
    return out;
  }

  /**
   * Computes log q and d(log q) / d(alpha) using log1m.
   */
  template <typename T, typename T_shape,
            bool partials_fvar>
  static inline Q_eval<T> eval_q_log1m(const T& alpha,
                                       const T& beta_y) {
    Q_eval<T> out;
    out.log_Q = log1m(gamma_p(alpha, beta_y));

    if (!std::isfinite(value_of_rec(out.log_Q))) {
      out.ok = false;
      return out;
    }

    if constexpr (is_autodiff_v<T_shape>) {
      if constexpr (partials_fvar) {
        T alpha_unit = alpha;
        alpha_unit.d_ = 1;
        T beta_unit = beta_y;
        beta_unit.d_ = 0;
        T log_Q_fvar = log1m(gamma_p(alpha_unit, beta_unit));
        out.dlogQ_dalpha = log_Q_fvar.d_;
      } else {
        out.dlogQ_dalpha = -grad_reg_lower_inc_gamma(alpha, beta_y) / exp(out.log_Q);
      }
    }

    out.ok = true;
    return out;
  }
}

template <typename T_y, typename T_shape, typename T_inv_scale>
inline return_type_t<T_y, T_shape, T_inv_scale> gamma_lccdf(const T_y& y,
                                                            const T_shape& alpha,
                                                            const T_inv_scale& beta) {
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
    const T_partials_return y_dbl = y_vec.val(n);
    if (y_dbl == 0.0) {
      continue;
    }
    if (y_dbl == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return alpha_dbl = alpha_vec.val(n);
    const T_partials_return beta_dbl = beta_vec.val(n);

    const T_partials_return beta_y = beta_dbl * y_dbl;
    if (beta_y == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const bool use_continued_fraction = beta_y > alpha_dbl + 1.0;
    internal::Q_eval<T_partials_return> result;
    if (use_continued_fraction) {
      result = internal::eval_q_cf<T_partials_return, T_shape,
                                any_fvar, partials_fvar>(alpha_dbl, beta_y);
    } else {
      result = internal::eval_q_log1m<T_partials_return, T_shape,
                                   partials_fvar>(alpha_dbl, beta_y);

      if (!result.ok && beta_y > 0.0) {
        // Fallback to continued fraction if log1m fails
        result = internal::eval_q_cf<T_partials_return, T_shape,
                                  any_fvar, partials_fvar>(alpha_dbl, beta_y);
      }
    }
    if (!result.ok) {
      return ops_partials.build(negative_infinity());
    }

    P += result.log_Q;

    if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_inv_scale>) {
      const T_partials_return log_y = log(y_dbl);
      const T_partials_return alpha_minus_one = fma(alpha_dbl, log_y, -log_y);

      const T_partials_return log_pdf
        = alpha_dbl * log(beta_dbl) - lgamma(alpha_dbl) + alpha_minus_one - beta_y;

      const T_partials_return hazard = exp(log_pdf - result.log_Q);  // f/Q

      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] -= hazard;
      }
      if constexpr (is_autodiff_v<T_inv_scale>) {
        partials<2>(ops_partials)[n] -= (y_dbl / beta_dbl) * hazard;
      }
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<1>(ops_partials)[n] += result.dlogQ_dalpha;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan

#endif
