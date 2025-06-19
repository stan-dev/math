#ifndef STAN_MATH_FWD_PROB_STUDENT_T_QF_HPP
#define STAN_MATH_FWD_PROB_STUDENT_T_QF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/fwd/fun/inv_inc_beta.hpp>
#include <stan/math/prim/prob/student_t_lpdf.hpp>

namespace stan {
namespace math {

template <typename T_p, typename T_nu, typename T_mu, typename T_sigma,
          typename FvarT = return_type_t<T_p, T_mu, T_sigma, T_nu>,
          require_all_stan_scalar_t<T_p, T_mu, T_sigma, T_nu>* = nullptr,
          require_fvar_t<FvarT>* = nullptr>
FvarT student_t_qf(const T_p& p, const T_nu& nu, const T_mu& mu,
                    const T_sigma& sigma) {
  static constexpr const char* function = "student_t_qf";
  using T_partials = partials_type_t<FvarT>;

  const auto& p_val = value_of(p);
  const auto& nu_val = value_of(nu);
  const auto& mu_val = value_of(mu);
  const auto& sigma_val = value_of(sigma);

  check_nonnegative(function, "Degrees of freedom parameter", nu_val);
  check_positive(function, "Scale parameter", sigma_val);
  check_bounded(function, "Probability parameter", p_val, 0.0, 1.0);

  if (p_val == 0.0) {
    return {NEGATIVE_INFTY, 0.0};
  } else if (p_val == 1.0) {
    return {INFTY, 0.0};
  } else if (p_val == 0.5) {
    return {mu_val, 0.0};
  }


  const T_partials p_val_flip = p_val < 0.5 ? p_val : 1.0 - p_val;
  const double p_sign = value_of_rec(p_val) < 0.5 ? -1.0 : 1.0;

  T_partials ibeta_arg = inv_inc_beta(0.5 * nu_val, 0.5, 2 * p_val_flip);
  T_partials rtn_val = mu_val + p_sign * sigma_val * sqrt(nu_val)
                * math::sqrt(-1.0 + 1.0 / ibeta_arg);

  FvarT rtn(rtn_val, 0.0);

  if constexpr (!is_constant<T_p>::value) {
    rtn.d_ += p.d_ * exp(-student_t_lpdf(rtn_val, nu_val, mu_val, sigma_val));
  }

  if constexpr (!is_constant<T_nu>::value) {
    const T_partials half_nu = nu_val / 2.0;
    Eigen::Matrix<T_partials, Eigen::Dynamic, 1> hyper_arg_a(3);
    hyper_arg_a << 0.5, half_nu, half_nu;
    Eigen::Matrix<T_partials, Eigen::Dynamic, 1> hyper_arg_b(2);
    hyper_arg_b << 1 + half_nu, 1 + half_nu;
    const T_partials hyper_arg = hypergeometric_pFq(hyper_arg_a, hyper_arg_b, ibeta_arg);
    const T_partials hyper2f1 = hypergeometric_2F1(1, (1 + nu_val) / 2, (2 + nu_val) / 2, ibeta_arg);
    const T_partials digamma_a1 = digamma(half_nu);
    const T_partials digamma_a2 = digamma((1 + nu_val) / 2);

    const T_partials arg_1 = (4 * hyper_arg * sqrt(1 - ibeta_arg)) / nu_val;
    const T_partials arg_2 = -2 * hyper2f1 * (-1 + ibeta_arg) * (log(ibeta_arg) - digamma_a1 + digamma_a2);

    const T_partials num1 = sigma_val * (-2 + (2 - arg_1 + arg_2) / ibeta_arg);
    const T_partials den1 = 4 * sqrt(nu_val) * sqrt(-1 + 1 / ibeta_arg);
    rtn.d_ += nu.d_ * p_sign * num1 / den1;
  }

  if constexpr (!is_constant<T_mu>::value) {
    rtn.d_ += mu.d_;
  }

  if constexpr (!is_constant<T_sigma>::value) {
    rtn.d_ += sigma.d_ * p_sign * sqrt(nu_val) * math::sqrt(-1.0 + 1.0 / ibeta_arg);
  }

  return rtn;

}
}
}

#endif
