#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_QF_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/inv_inc_beta.hpp>
#include <stan/math/prim/fun/max_size.hpp>

namespace stan {
namespace math {

/**
 * The quantile function of the Student's t-distribution.
 *
 * @tparam T_p type of the probability parameter
 * @tparam T_nu type of the degrees of freedom parameter
 * @tparam T_mu type of the location parameter
 * @tparam T_sigma type of the scale parameter
 * @param p Probability in the range [0, 1].
 * @param nu Degrees of freedom, must be non-negative.
 * @param mu Location parameter.
 * @param sigma Scale parameter, must be positive.
 * @return Quantile function value.
 * @throw std::domain_error if `nu` is negative or `sigma` is not positive,
 * or if `p` is not in [0, 1].
 */
template <typename T_p, typename T_nu, typename T_mu, typename T_sigma,
          require_all_stan_scalar_t<T_p, T_nu, T_mu, T_sigma>* = nullptr,
          require_all_arithmetic_t<T_p, T_nu, T_mu, T_sigma>* = nullptr>
inline double student_t_qf(const T_p& p, const T_nu& nu, const T_mu& mu,
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
  const auto ibeta_arg = inv_inc_beta(0.5 * nu, 0.5, 2 * p_val_flip);

  return mu + p_sign * sigma * sqrt(nu) * sqrt(-1.0 + 1.0 / ibeta_arg);
}

/**
 * A vectorized version of the Student's t quantile function that accepts
 * std::vectors, Eigen Matrix/Array objects, or expressions, and containers of
 * these.
 *
 * @tparam T_p type of the probability parameter
 * @tparam T_nu type of the degrees of freedom parameter
 * @tparam T_mu type of the location parameter
 * @tparam T_sigma type of the scale parameter
 * @tparam T_container type of the container to hold results
 * @param p Probability in the range [0, 1].
 * @param nu Degrees of freedom, must be non-negative.
 * @param mu Location parameter.
 * @param sigma Scale parameter, must be positive.
 * @return Container with quantile function values for each input.
 */
template <typename T_p, typename T_nu, typename T_mu, typename T_sigma,
          require_any_vector_t<T_p, T_nu, T_mu, T_sigma>* = nullptr>
inline auto student_t_qf(const T_p& p, const T_nu& nu, const T_mu& mu,
                         const T_sigma& sigma) {
  using T_container = common_container_t<T_p, T_nu, T_mu, T_sigma>;
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

}  // namespace math
}  // namespace stan

#endif
