#ifndef STAN_MATH_FWD_PROB_STUDENT_T_QF_HPP
#define STAN_MATH_FWD_PROB_STUDENT_T_QF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/hypergeometric_2F1.hpp>
#include <stan/math/fwd/fun/hypergeometric_pFq.hpp>
#include <stan/math/fwd/fun/inv_inc_beta.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/fwd/fun/value_of_rec.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/student_t_lpdf.hpp>

namespace stan {
namespace math {

template <typename T_p, typename T_nu, typename T_mu, typename T_sigma,
          require_all_stan_scalar_t<T_p, T_mu, T_sigma, T_nu>* = nullptr,
          require_any_fvar_t<T_p, T_nu, T_mu, T_sigma>* = nullptr>
inline auto student_t_qf(const T_p& p, const T_nu& nu, const T_mu& mu,
                          const T_sigma& sigma) {
  static constexpr const char* function = "student_t_qf";
  using FvarT = return_type_t<T_p, T_mu, T_sigma, T_nu>;
  using T_partials = partials_type_t<FvarT>;

  auto p_val = value_of(p);
  auto nu_val = value_of(nu);
  auto mu_val = value_of(mu);
  auto sigma_val = value_of(sigma);

  check_nonnegative(function, "Degrees of freedom parameter", nu_val);
  check_positive(function, "Scale parameter", sigma_val);
  check_bounded(function, "Probability parameter", p_val, 0.0, 1.0);

  if (unlikely(p_val == 0.0)) {
    return FvarT{NEGATIVE_INFTY, 0.0};
  } else if (unlikely(p_val == 1.0)) {
    return FvarT{INFTY, 0.0};
  } else if (unlikely(p_val == 0.5)) {
    return FvarT{mu_val, 0.0};
  }

  const auto p_val_flip = p_val < 0.5 ? p_val : 1.0 - p_val;
  const double p_sign = value_of_rec(p_val) < 0.5 ? -1.0 : 1.0;
  auto sqrt_nu_val = sqrt(nu_val);
  auto ibeta_arg = inv_inc_beta(0.5 * nu_val, 0.5, 2.0 * p_val_flip);
  auto rtn_val = mu_val
                 + p_sign * sigma_val * sqrt_nu_val
                       * sqrt(-1.0 + 1.0 / ibeta_arg);

  FvarT rtn(rtn_val, 0.0);

  if constexpr (is_autodiff_v<T_p>) {
    rtn.d_ += p.d_ * exp(-student_t_lpdf(rtn_val, nu_val, mu_val, sigma_val));
  }

  if constexpr (is_autodiff_v<T_nu>) {
    const auto half_nu = nu_val / 2.0;
    Eigen::Matrix<T_partials, -1, 1> hyper_arg_a{{0.5, half_nu, half_nu}};
    Eigen::Matrix<T_partials, -1, 1> hyper_arg_b{{1.0 + half_nu, 1.0 + half_nu}};
    const auto hyper_arg
        = hypergeometric_pFq(hyper_arg_a, hyper_arg_b, ibeta_arg);
    const auto hyper2f1
        = hypergeometric_2F1(1.0, (1.0 + nu_val) / 2.0, (2.0 + nu_val) / 2.0, ibeta_arg);
    const auto digamma_a1 = digamma(half_nu);
    const auto digamma_a2 = digamma((1.0 + nu_val) / 2.0);

    const auto arg_1 = (4.0 * hyper_arg * sqrt(1.0 - ibeta_arg)) / nu_val;
    const auto arg_2 = -2.0 * hyper2f1 * (-1.0 + ibeta_arg)
                             * (log(ibeta_arg) - digamma_a1 + digamma_a2);

    const auto num1 = sigma_val * (-2.0 + (2.0 - arg_1 + arg_2) / ibeta_arg);
    const auto den1 = 4.0 * sqrt_nu_val * sqrt(-1.0 + 1.0 / ibeta_arg);
    rtn.d_ += nu.d_ * p_sign * num1 / den1;
  }

  if constexpr (is_autodiff_v<T_mu>) {
    rtn.d_ += mu.d_;
  }

  if constexpr (is_autodiff_v<T_sigma>) {
    rtn.d_ += sigma.d_ * p_sign * sqrt_nu_val
              * sqrt(-1.0 + 1.0 / ibeta_arg);
  }

  return rtn;
}
}  // namespace math
}  // namespace stan

#endif
