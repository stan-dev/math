#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_LCDF_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_dof, T_loc, T_scale> student_t_lcdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_dof, T_loc, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  using std::exp;
  using std::log;
  using std::pow;
  static constexpr const char* function = "student_t_lcdf";
  T_y_ref y_ref = y;
  T_nu_ref nu_ref = nu;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_positive_finite(function, "Degrees of freedom parameter", nu_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  if (size_zero(y, nu, mu, sigma)) {
    return 0;
  }

  T_partials_return P(0.0);
  auto ops_partials
      = make_partials_propagator(y_ref, nu_ref, mu_ref, sigma_ref);
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_nu_ref> nu_vec(nu_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t N = max_size(y, nu, mu, sigma);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) == NEGATIVE_INFTY) {
      return ops_partials.build(negative_infinity());
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_vec.val(n) == INFTY) {
      continue;
    }

    const auto sigma_inv = 1.0 / sigma_vec.val(n);
    const auto t = (y_vec.val(n) - mu_vec.val(n)) * sigma_inv;
    const auto nu_dbl = nu_vec.val(n);
    const auto q = nu_dbl / (t * t);
    const auto r = 1.0 / (1.0 + q);
    // 1 - r = 1 / (1 + 1 / q) without the subtraction, which rounds to 0
    // once q is below eps; equal to 1 at t = 0
    const auto one_m_r = 1.0 / (1.0 + t * t / nu_dbl);
    const auto J = 2 * r * r * q / t;
    const auto betaNuHalf = beta(0.5, 0.5 * nu_dbl);
    auto zJacobian = t > 0 ? -0.5 : 0.5;

    if (q < 2) {
      auto z = inc_beta(0.5 * nu_dbl, 0.5, one_m_r);
      const auto Pn = t > 0 ? 1.0 - 0.5 * z : 0.5 * z;
      const auto d_ibeta
          = pow(r, -0.5) * pow(one_m_r, 0.5 * nu_dbl - 1) / betaNuHalf;

      P += log(Pn);

      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n]
            += -zJacobian * d_ibeta * J * sigma_inv / Pn;
      }

      if constexpr (is_autodiff_v<T_dof>) {
        // nu / 2 is the first shape of inc_beta here and 1 / 2 the second
        const auto g1 = grad_reg_inc_beta<GRAD_RIB::G1>(
            0.5 * nu_dbl, 0.5, one_m_r, digamma(0.5 * nu_dbl), 0.0,
            digamma(0.5 + 0.5 * nu_dbl), 0.0);

        partials<1>(ops_partials)[n]
            += zJacobian * (d_ibeta * (r / t) * (r / t) + 0.5 * g1) / Pn;
      }

      if constexpr (is_autodiff_v<T_loc>) {
        partials<2>(ops_partials)[n]
            += zJacobian * d_ibeta * J * sigma_inv / Pn;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<3>(ops_partials)[n]
            += zJacobian * d_ibeta * J * sigma_inv * t / Pn;
      }

    } else {
      // Same quantity as in the other branch, by the symmetry relation
      // 1 - I_r(1/2, nu/2) = I_{1-r}(nu/2, 1/2). Forming it as the
      // complement loses every digit once I_r(1/2, nu/2) rounds to 1,
      // which happens for every z below eps, that is for every tail
      // probability below about e^-37.
      auto z = inc_beta(0.5 * nu_dbl, 0.5, one_m_r);
      zJacobian *= -1;

      const auto Pn = t > 0 ? 1.0 - 0.5 * z : 0.5 * z;

      auto d_ibeta = pow(one_m_r, 0.5 * nu_dbl - 1) * pow(r, -0.5) / betaNuHalf;

      P += log(Pn);

      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n]
            += zJacobian * d_ibeta * J * sigma_inv / Pn;
      }

      if constexpr (is_autodiff_v<T_dof>) {
        // On the reflected arguments nu / 2 is the first shape, and the
        // partial is d z / d (nu / 2) directly, which is the negative of the
        // second-shape partial of the unreflected call.
        const auto g1 = grad_reg_inc_beta<GRAD_RIB::G1>(
            0.5 * nu_dbl, 0.5, one_m_r, digamma(0.5 * nu_dbl), 0.0,
            digamma(0.5 + 0.5 * nu_dbl), 0.0);

        partials<1>(ops_partials)[n]
            += zJacobian * (-d_ibeta * (r / t) * (r / t) - 0.5 * g1) / Pn;
      }

      if constexpr (is_autodiff_v<T_loc>) {
        partials<2>(ops_partials)[n]
            += -zJacobian * d_ibeta * J * sigma_inv / Pn;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<3>(ops_partials)[n]
            += -zJacobian * d_ibeta * J * sigma_inv * t / Pn;
      }
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
