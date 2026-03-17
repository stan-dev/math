#ifndef STAN_MATH_REV_PROB_STUDENT_T_QF_HPP
#define STAN_MATH_REV_PROB_STUDENT_T_QF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/hypergeometric_pFq.hpp>
#include <stan/math/rev/fun/hypergeometric_2F1.hpp>
#include <stan/math/rev/fun/inv.hpp>
#include <stan/math/rev/fun/inv_inc_beta.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/prob/student_t_lpdf.hpp>

namespace stan {
namespace math {

template <typename T_p, typename T_nu, typename T_mu, typename T_sigma,
          require_all_stan_scalar_t<T_p, T_mu, T_sigma, T_nu>* = nullptr,
          require_any_var_t<T_p, T_mu, T_sigma, T_nu>* = nullptr>
inline var student_t_qf(const T_p& p, const T_nu& nu, const T_mu& mu,
                        const T_sigma& sigma) {
  static constexpr const char* function = "student_t_qf";
  const double p_val = value_of(p);
  const double nu_val = value_of(nu);
  const double mu_val = value_of(mu);
  const double sigma_val = value_of(sigma);
  check_nonnegative(function, "Degrees of freedom parameter", nu_val);
  check_positive(function, "Scale parameter", sigma_val);
  check_bounded(function, "Probability parameter", p_val, 0.0, 1.0);
  if (unlikely(p == 0.0)) {
    return var{NEGATIVE_INFTY};
  } else if (unlikely(p == 1.0)) {
    return var{INFTY};
  } else if (unlikely(p == 0.5)) {
    return var{mu_val};
  } else {
    const double p_val_flip = p_val < 0.5 ? p_val : 1.0 - p_val;
    const double p_sign = p_val < 0.5 ? -1.0 : 1.0;
    const double ibeta_arg = inv_inc_beta(0.5 * nu_val, 0.5, 2 * p_val_flip);
    const double sqrt_inv_ibeta_m1 = sqrt(inv(ibeta_arg) - 1.0);
    double ret_val = mu_val
              + p_sign * sigma_val * sqrt(nu_val)
                    * sqrt(-1.0 + 1.0 / ibeta_arg);
    return make_callback_var(ret_val, [p, mu, sigma, nu,
                                      ibeta_arg](auto& vi) mutable {
      const double p_val = value_of(p);
      const double mu_val = value_of(mu);
      const double sigma_val = value_of(sigma);
      const double nu_val = value_of(nu);
      const double p_val_flip = p_val < 0.5 ? p_val : 1.0 - p_val;
      const double p_sign = p_val < 0.5 ? -1.0 : 1.0;
      const double sqrt_inv_ibeta_m1 = sqrt(inv(ibeta_arg) - 1.0);
      if constexpr (is_autodiff_v<T_p>) {
        if (likely(p.val() > 0.0 && p.val() < 1.0)) {
          p.adj() += vi.adj()
                    * exp(-student_t_lpdf(vi.val(), nu_val, mu_val, sigma_val));
        }
      }
      if constexpr (is_autodiff_v<T_nu>) {
        const double half_nu = nu_val / 2.0;
        const double nu_p1_div_2 = (1.0 + nu_val) / 2.0;
        Eigen::VectorXd hyper_arg_a{{0.5, half_nu, half_nu}};
        Eigen::VectorXd hyper_arg_b{{1.0 + half_nu, 1.0 + half_nu}};
        const double hyper_arg
            = hypergeometric_pFq(hyper_arg_a, hyper_arg_b, ibeta_arg);
        const double hyper2f1 = hypergeometric_2F1(
            1.0, nu_p1_div_2, (2.0 + nu_val) / 2.0, ibeta_arg);
        const double digamma_a1 = digamma(half_nu);
        const double digamma_a2 = digamma(nu_p1_div_2);

        const double arg_1 = (4.0 * hyper_arg * sqrt(1.0 - ibeta_arg)) / nu_val;
        const double arg_2 = -2.0 * hyper2f1 * (-1.0 + ibeta_arg)
                            * (log(ibeta_arg) - digamma_a1 + digamma_a2);

        const double num1 = sigma_val * (-2.0 + (2.0 - arg_1 + arg_2) / ibeta_arg);
        const double den1 = 4.0 * sqrt(nu_val) * sqrt_inv_ibeta_m1;
        nu.adj() += vi.adj() * p_sign * num1 / den1;
      }
      if constexpr (is_autodiff_v<T_mu>) {
        if (p_val > 0.0 && p_val < 1.0) {
          mu.adj() += vi.adj();
        }
      }
      if constexpr (is_autodiff_v<T_sigma>) {
        if (p_val > 0.0 && p_val < 1.0) {
          sigma.adj() += vi.adj() * p_sign * sqrt(nu_val) * sqrt_inv_ibeta_m1;
        }
      }
    });
  }
}
}  // namespace math
}  // namespace stan

#endif
