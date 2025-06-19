#ifndef STAN_MATH_REV_PROB_STUDENT_T_QF_HPP
#define STAN_MATH_REV_PROB_STUDENT_T_QF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <stan/math/rev/fun/inv_inc_beta.hpp>

namespace stan {
namespace math {

template <typename T_p, typename T_mu, typename T_sigma, typename T_nu,
          require_all_stan_scalar_t<T_p, T_mu, T_sigma, T_nu>* = nullptr,
          require_any_var_t<T_p, T_mu, T_sigma, T_nu>* = nullptr>
var student_t_qf(const T_p& p, const T_mu& mu,
                    const T_sigma& sigma, const T_nu& nu) {
  static constexpr const char* function = "student_t_qf";
  const double p_val = value_of(p);
  const double mu_val = value_of(mu);
  const double sigma_val = value_of(sigma);
  const double nu_val = value_of(nu);

  check_nonnegative(function, "Degrees of freedom parameter", nu_val);
  check_positive(function, "Scale parameter", sigma_val);
  check_bounded(function, "Probability parameter", p_val, 0.0, 1.0);

  double ret_val(0);
  double ibeta_arg(0);

  if (p == 0.0) {
    ret_val = NEGATIVE_INFTY;
  } else if (p == 1.0) {
    ret_val = INFTY;
  } else if (p == 0.5) {
    ret_val = mu_val;
  } else {
    const double p_val_flip = p_val < 0.5 ? p_val : 1.0 - p_val;
    const double p_sign = p_val < 0.5 ? -1.0 : 1.0;
    ibeta_arg = inv_inc_beta(0.5 * nu_val, 0.5, 2 * p_val_flip);

    ret_val = mu_val + p_sign * sigma_val * math::sqrt(nu_val)
                  * math::sqrt(-1.0 + 1.0 / ibeta_arg);
  }

  return make_callback_var(ret_val, [p, mu, sigma, nu, ibeta_arg, nu_val, ret_val](auto& vi) mutable {
    const double p_val = value_of(p);
    if constexpr (!is_constant<T_p>::value) {
      if (p.val() > 0.0 && p.val() < 1.0) {
        p.adj() += vi.adj() * exp(log(p_val) - student_t_lpdf(ret_val, mu.val(), sigma.val(), nu_val));
      }
    }
    if constexpr (!is_constant<T_mu>::value) {
      if (p_val > 0.0 && p_val < 1.0) {
        mu.adj() += vi.adj();
      }
    }
    if constexpr (!is_constant<T_sigma>::value) {
      if (p_val > 0.0 && p_val < 1.0) {
        sigma.adj() += vi.adj();
      }
    }
    if constexpr (!is_constant<T_nu>::value) {
      const double sigma_val = value_of(sigma);
      const double p_sign = p_val < 0.5 ? -1.0 : 1.0;
      const double half_nu = nu_val / 2.0;
      const double hyper_arg = hypergeometric_3F2({0.5, half_nu, half_nu}, {1 + half_nu, 1 + half_nu}, ibeta_arg);
      const double hyper2f1 = hypergeometric_2F1(1, (1 + nu_val) / 2, (2 + nu_val) / 2, ibeta_arg);
      const double digamma_a1 = digamma(half_nu);
      const double digamma_a2 = digamma((1 + nu_val) / 2);

      const double arg_1 = (4 * hyper_arg * sqrt(1 - ibeta_arg)) / nu_val;
      const double arg_2 = -2 * hyper2f1 * (-1 + ibeta_arg) * (log(ibeta_arg) - digamma_a1 + digamma_a2);

      const double num1 = sigma_val * (-2 + (2 - arg_1 + arg_2) / ibeta_arg);
      const double den1 = 4 * sqrt(nu_val) * sqrt(-1 + 1 / ibeta_arg);
      nu.adj() += vi.adj() * p_sign * num1 / den1;
    }
  });

}
}
}

#endif
