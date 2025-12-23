#ifndef STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fma.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/grad_reg_lower_inc_gamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/fwd/fun/lgamma.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

template <typename T_a>
inline auto log_q_gamma_cf(const T_a& a, const double x,
                           int max_steps = 250, double precision = 1e-16) {
  using std::fabs;
  using stan::math::lgamma;
  using stan::math::log;
  using stan::math::value_of;

  const auto log_prefactor = a * log(x) - x - lgamma(a);

  auto b = x + 1.0 - a;
  auto C = (fabs(value_of(b)) >= EPSILON) ? b : T_a(EPSILON);
  auto D = T_a(0.0);
  auto f = C;

  for (int i = 1; i <= max_steps; ++i) {
    auto an = -i * (i - a);
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

    const double delta_m1 = fabs(value_of(delta) - 1.0);
    if (delta_m1 < precision) {
      if constexpr (stan::is_fvar<std::decay_t<T_a>>::value) {
        if (fabs(value_of(delta.d_)) < precision) {
          break;
        }
      } else {
        break;
      }
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

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    const T_partials_return y_dbl = value_of(y_vec.val(n));
    if (y_dbl == 0.0) {
      continue;
    }
    if (y_dbl == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return alpha_dbl = value_of(alpha_vec.val(n));
    const T_partials_return beta_dbl = value_of(beta_vec.val(n));

    const T_partials_return beta_y = beta_dbl * y_dbl;
    if (beta_y == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    bool use_cf = beta_y > alpha_dbl + 1.0;
    T_partials_return log_Qn;
    [[maybe_unused]] T_partials_return dlogQ_dalpha = 0.0;
    // Extract double values for continued fraction - we handle y/beta derivatives via hazard
    const double beta_y_dbl = value_of(value_of(beta_y));
    const double alpha_dbl_val = value_of(value_of(alpha_dbl));

    if (use_cf) {
      if constexpr (is_autodiff_v<T_shape>) {
        stan::math::fvar<double> a_f(alpha_dbl_val, 1.0);
        const stan::math::fvar<double> logq_f
            = internal::log_q_gamma_cf(a_f, beta_y_dbl);
        log_Qn = logq_f.val_;
        dlogQ_dalpha = logq_f.d_;
      } else {
        log_Qn = internal::log_q_gamma_cf(alpha_dbl_val, beta_y_dbl);
      }
    } else {
      const T_partials_return Pn = gamma_p(alpha_dbl, beta_y);
      log_Qn = log1m(Pn);

      if (!std::isfinite(value_of(value_of(log_Qn)))) {
        use_cf = beta_y > 0.0;
        if (use_cf) {
          if constexpr (is_autodiff_v<T_shape>) {
            stan::math::fvar<double> a_f(alpha_dbl_val, 1.0);
            const stan::math::fvar<double> logq_f
                = internal::log_q_gamma_cf(a_f, beta_y_dbl);
            log_Qn = logq_f.val_;
            dlogQ_dalpha = logq_f.d_;
          } else {
            log_Qn = internal::log_q_gamma_cf(alpha_dbl_val, beta_y_dbl);
          }
        }
      }

      if constexpr (is_autodiff_v<T_shape>) {
        if (!use_cf) {
          const T_partials_return Qn = exp(log_Qn);
          if (Qn > 0.0) {
            dlogQ_dalpha = -grad_reg_lower_inc_gamma(alpha_dbl, beta_y) / Qn;
          } else {
            stan::math::fvar<double> a_f(alpha_dbl_val, 1.0);
            const stan::math::fvar<double> logq_f
                = internal::log_q_gamma_cf(a_f, beta_y_dbl);
            log_Qn = logq_f.val_;
            dlogQ_dalpha = logq_f.d_;
            use_cf = true;
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

      const T_partials_return log_pdf = alpha_dbl * log_beta - lgamma_alpha
                                        + alpha_minus_one - beta_y;

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
