#ifndef STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_dof, typename T_scale>
return_type_t<T_y, T_dof, T_scale> scaled_inv_chi_square_lccdf(
    const T_y& y, const T_dof& nu, const T_scale& s) {
  using T_partials_return = partials_return_t<T_y, T_dof, T_scale>;
  using std::exp;
  using std::log;
  using std::pow;
  static const char* function = "scaled_inv_chi_square_lccdf";
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_positive_finite(function, "Scale parameter", s);
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu, "Scale parameter",
                         s);

  if (size_zero(y, nu, s)) {
    return 0;
  }

  T_partials_return P(0.0);
  operands_and_partials<T_y, T_dof, T_scale> ops_partials(y, nu, s);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_dof> nu_vec(nu);
  scalar_seq_view<T_scale> s_vec(s);
  size_t size_y = stan::math::size(y);
  size_t size_nu = stan::math::size(nu);
  size_t N = max_size(y, nu, s);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size_y; i++) {
    const T_partials_return y_dbl = value_of(y_vec[i]);
    if (y_dbl == 0) {
      return ops_partials.build(0.0);
    }
    if (y_dbl == INFTY) {
      return ops_partials.build(NEGATIVE_INFTY);
    }
  }

  VectorBuilder<true, T_partials_return, T_y> inv_y(size_y);
  for (size_t i = 0; i < size_y; i++) {
    inv_y[i] = inv(value_of(y_vec[i]));
  }

  VectorBuilder<true, T_partials_return, T_dof> tgamma_vec(size_nu);
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_vec(size_nu);
  for (size_t i = 0; i < size_nu; i++) {
    const T_partials_return half_nu_dbl = 0.5 * value_of(nu_vec[i]);
    tgamma_vec[i] = tgamma(half_nu_dbl);
    if (!is_constant_all<T_dof>::value) {
      digamma_vec[i] = digamma(half_nu_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return half_nu_dbl = 0.5 * value_of(nu_vec[n]);
    const T_partials_return s_dbl = value_of(s_vec[n]);
    const T_partials_return s2_over_y = square(s_dbl) * inv_y[n];
    const T_partials_return half_nu_s2_over_y = half_nu_dbl * s2_over_y;
    const T_partials_return Pn = gamma_p(half_nu_dbl, half_nu_s2_over_y);
    const T_partials_return gamma_p_deriv
        = exp(-half_nu_s2_over_y) * pow(half_nu_s2_over_y, half_nu_dbl - 1)
          / tgamma_vec[n];

    P += log(Pn);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          -= half_nu_s2_over_y * inv_y[n] * gamma_p_deriv / Pn;
    }
    if (!is_constant_all<T_dof>::value) {
      ops_partials.edge2_.partials_[n]
          -= 0.5
             * (grad_reg_inc_gamma(half_nu_dbl, half_nu_s2_over_y,
                                   tgamma_vec[n], digamma_vec[n])
                - s2_over_y * gamma_p_deriv)
             / Pn;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += 2.0 * half_nu_dbl * s_dbl * inv_y[n] * gamma_p_deriv / Pn;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
