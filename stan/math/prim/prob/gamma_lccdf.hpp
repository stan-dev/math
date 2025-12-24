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
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/log_gamma_q_dgamma.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

/**
 * Compute log(Q(a,x)) using continued fraction expansion for upper incomplete
 * gamma function. When used with fvar types, automatically computes
 * derivatives.
 *
 * @tparam T_a Type of shape parameter a (double or fvar types)
 * @param a Shape parameter
 * @param x Value at which to evaluate
 * @param max_steps Maximum number of continued fraction iterations
 * @param precision Convergence threshold
 * @return log(Q(a,x)) with same type as T_a
 */
template <typename T_a, typename T_x>
inline auto log_q_gamma_cf(const T_a& a, const T_x& x, int max_steps = 250,
                           double precision = 1e-16) {
  using stan::math::lgamma;
  using stan::math::log;
  using stan::math::value_of;
  using std::fabs;
  using T_return = return_type_t<T_a, T_x>;

  const T_return a_ret = a;
  const T_return x_ret = x;
  const auto log_prefactor = a_ret * log(x_ret) - x_ret - lgamma(a_ret);

  auto b = x_ret + 1.0 - a_ret;
  auto C = (fabs(value_of(b)) >= EPSILON) ? b : T_return(EPSILON);
  auto D = T_return(0.0);
  auto f = C;

  for (int i = 1; i <= max_steps; ++i) {
    auto an = -i * (i - a_ret);
    b += 2.0;

    D = b + an * D;
    if (fabs(value_of(D)) < EPSILON) {
      D = T_a(EPSILON);
    }
    C = b + an / C;
    if (fabs(value_of(C)) < EPSILON) {
      C = T_a(EPSILON);
    }

    D = 1.0 / D;
    auto delta = C * D;
    f *= delta;

    const double delta_m1 = value_of(fabs(value_of(delta) - 1.0));
    if (delta_m1 < precision) {
      break;
    }
  }

  return log_prefactor - log(f);
}

}  // namespace internal

template <typename T_y, typename T_shape, typename T_inv_scale>
return_type_t<T_y, T_shape, T_inv_scale> gamma_lccdf(const T_y& y,
                                                     const T_shape& alpha,
                                                     const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using std::exp;
  using std::log;
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

  constexpr bool need_y_beta_deriv = !is_constant_all<T_y, T_inv_scale>::value;
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

    bool use_cf = beta_y > alpha_dbl + 1.0;
    T_partials_return log_Qn;
    [[maybe_unused]] T_partials_return dlogQ_dalpha = 0.0;
    // Extract double values for the double-only continued fraction path.
    [[maybe_unused]] const double beta_y_dbl = value_of(value_of(beta_y));
    [[maybe_unused]] const double alpha_dbl_val = value_of(value_of(alpha_dbl));

    if (use_cf) {
      if constexpr (!any_fvar && is_autodiff_v<T_shape>) {
        // var-only: use analytical gradient with double inputs
        auto log_q_result = log_gamma_q_dgamma(alpha_dbl_val, beta_y_dbl);
        log_Qn = log_q_result.log_q;
        dlogQ_dalpha = log_q_result.dlog_q_da;
      } else {
        log_Qn = internal::log_q_gamma_cf(alpha_dbl, beta_y);
        if constexpr (is_autodiff_v<T_shape>) {
          if constexpr (partials_fvar) {
            auto alpha_unit = alpha_dbl;
            alpha_unit.d_ = 1;
            auto beta_unit = beta_y;
            beta_unit.d_ = 0;
            auto log_Qn_fvar = internal::log_q_gamma_cf(alpha_unit, beta_unit);
            dlogQ_dalpha = log_Qn_fvar.d_;
          } else {
            const T_partials_return Qn = exp(log_Qn);
            dlogQ_dalpha
                = grad_reg_inc_gamma(alpha_dbl, beta_y, tgamma(alpha_dbl),
                                     digamma(alpha_dbl))
                  / Qn;
          }
        }
      }
    } else {
      const T_partials_return Pn = gamma_p(alpha_dbl, beta_y);
      log_Qn = log1m(Pn);

      if (!std::isfinite(value_of(value_of(log_Qn)))) {
        use_cf = beta_y > 0.0;
        if (use_cf) {
          // Fallback to continued fraction if log1m fails
          if constexpr (!any_fvar && is_autodiff_v<T_shape>) {
            auto log_q_result = log_gamma_q_dgamma(alpha_dbl_val, beta_y_dbl);
            log_Qn = log_q_result.log_q;
            dlogQ_dalpha = log_q_result.dlog_q_da;
          } else {
            log_Qn = internal::log_q_gamma_cf(alpha_dbl, beta_y);
            if constexpr (is_autodiff_v<T_shape>) {
              if constexpr (partials_fvar) {
                auto alpha_unit = alpha_dbl;
                alpha_unit.d_ = 1;
                auto beta_unit = beta_y;
                beta_unit.d_ = 0;
                auto log_Qn_fvar
                    = internal::log_q_gamma_cf(alpha_unit, beta_unit);
                dlogQ_dalpha = log_Qn_fvar.d_;
              } else {
                const T_partials_return Qn = exp(log_Qn);
                dlogQ_dalpha
                    = grad_reg_inc_gamma(alpha_dbl, beta_y, tgamma(alpha_dbl),
                                         digamma(alpha_dbl))
                      / Qn;
              }
            }
          }
        }
      }

      if constexpr (is_autodiff_v<T_shape>) {
        if (!use_cf) {
          if constexpr (partials_fvar) {
            auto alpha_unit = alpha_dbl;
            alpha_unit.d_ = 1;
            auto beta_unit = beta_y;
            beta_unit.d_ = 0;
            auto log_Qn_fvar = log1m(gamma_p(alpha_unit, beta_unit));
            dlogQ_dalpha = log_Qn_fvar.d_;
          } else {
            const T_partials_return Qn = exp(log_Qn);
            if (Qn > 0.0) {
              dlogQ_dalpha = -grad_reg_lower_inc_gamma(alpha_dbl, beta_y) / Qn;
            } else {
              // Fallback to continued fraction if Q rounds to zero
              if constexpr (!any_fvar) {
                auto log_q_result
                    = log_gamma_q_dgamma(alpha_dbl_val, beta_y_dbl);
                log_Qn = log_q_result.log_q;
                dlogQ_dalpha = log_q_result.dlog_q_da;
              } else {
                log_Qn = internal::log_q_gamma_cf(alpha_dbl, beta_y);
                const T_partials_return Qn_cf = exp(log_Qn);
                dlogQ_dalpha
                    = grad_reg_inc_gamma(alpha_dbl, beta_y, tgamma(alpha_dbl),
                                         digamma(alpha_dbl))
                      / Qn_cf;
              }
              use_cf = true;
            }
          }
        }
      }
    }
    if (!std::isfinite(value_of(value_of(log_Qn)))) {
      return ops_partials.build(negative_infinity());
    }
    P += log_Qn;

    if constexpr (need_y_beta_deriv) {
      const T_partials_return log_y = log(y_dbl);
      const T_partials_return log_beta = log(beta_dbl);
      const T_partials_return lgamma_alpha = lgamma(alpha_dbl);
      const T_partials_return alpha_minus_one = fma(alpha_dbl, log_y, -log_y);

      const T_partials_return log_pdf
          = alpha_dbl * log_beta - lgamma_alpha + alpha_minus_one - beta_y;

      const T_partials_return hazard = exp(log_pdf - log_Qn);  // f/Q

      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] -= hazard;
      }
      if constexpr (is_autodiff_v<T_inv_scale>) {
        partials<2>(ops_partials)[n] -= (y_dbl / beta_dbl) * hazard;
      }
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<1>(ops_partials)[n] += dlogQ_dalpha;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan

#endif
