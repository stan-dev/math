#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_QF_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/inv_inc_beta.hpp>
#include <stan/math/prim/fun/max_size.hpp>

namespace stan {
namespace math {

template <typename T_p, typename T_nu, typename T_mu, typename T_sigma,
          require_all_stan_scalar_t<T_p, T_nu, T_mu, T_sigma>* = nullptr,
          require_all_arithmetic_t<T_p, T_nu, T_mu, T_sigma>* = nullptr>
double student_t_qf(const T_p& p, const T_nu& nu, const T_mu& mu,
                    const T_sigma& sigma) {
  static constexpr const char* function = "student_t_qf";
  check_nonnegative(function, "Degrees of freedom parameter", nu);
  check_positive(function, "Scale parameter", sigma);
  check_bounded(function, "Probability parameter", p, 0.0, 1.0);

  if (p == 0.0) {
    return NEGATIVE_INFTY;
  } else if (p == 1.0) {
    return INFTY;
  } else if (p == 0.5) {
    return mu;
  }

  const double p_val_flip = p < 0.5 ? p : 1.0 - p;
  const double p_sign = p < 0.5 ? -1.0 : 1.0;
  const auto& ibeta_arg = inv_inc_beta(0.5 * nu, 0.5, 2 * p_val_flip);

  return mu + p_sign * sigma * sqrt(nu) * sqrt(-1.0 + 1.0 / ibeta_arg);
}

template <typename T_p, typename T_nu, typename T_mu, typename T_sigma,
          typename T_container = common_container_t<T_p, T_nu, T_mu, T_sigma>,
          require_any_vector_t<T_p, T_nu, T_mu, T_sigma>* = nullptr,
          require_not_t<std::is_void<T_container>>* = nullptr>
auto student_t_qf(const T_p& p, const T_nu& nu, const T_mu& mu,
                  const T_sigma& sigma) {
  static constexpr const char* function = "student_t_qf";
  const size_t max_size_all = max_size(p, nu, mu, sigma);
  T_container result(max_size_all);

  ref_type_t<T_p> p_ref = p;
  ref_type_t<T_nu> nu_ref = nu;
  ref_type_t<T_mu> mu_ref = mu;
  ref_type_t<T_sigma> sigma_ref = sigma;

  scalar_seq_view<ref_type_t<T_p>> p_vec(p_ref);
  scalar_seq_view<ref_type_t<T_nu>> nu_vec(nu_ref);
  scalar_seq_view<ref_type_t<T_mu>> mu_vec(mu_ref);
  scalar_seq_view<ref_type_t<T_sigma>> sigma_vec(sigma_ref);

  for (size_t i = 0; i < max_size_all; ++i) {
    result[i] = student_t_qf(p_vec[i], nu_vec[i], mu_vec[i], sigma_vec[i]);
  }

  return result;
}

}
}

#endif
