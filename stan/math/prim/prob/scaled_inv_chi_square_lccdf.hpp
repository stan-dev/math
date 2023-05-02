#ifndef STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
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
  using T_y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_s_ref = ref_type_t<T_scale>;
  static const char* function = "scaled_inv_chi_square_lccdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu, "Scale parameter",
                         s);
  T_y_ref y_ref = y;
  T_nu_ref nu_ref = nu;
  T_s_ref s_ref = s;
  check_nonnegative(function, "Random variable", y_ref);
  check_positive_finite(function, "Degrees of freedom parameter", nu_ref);
  check_positive_finite(function, "Scale parameter", s_ref);

  if (size_zero(y, nu, s)) {
    return 0;
  }

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(y_ref, nu_ref, s_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_nu_ref> nu_vec(nu_ref);
  scalar_seq_view<T_s_ref> s_vec(s_ref);
  size_t N = max_size(y, nu, s);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) == 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      gamma_vec(math::size(nu));
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_vec(math::size(nu));

  if (!is_constant_all<T_dof>::value) {
    for (size_t i = 0; i < stan::math::size(nu); i++) {
      const T_partials_return half_nu_dbl = 0.5 * nu_vec.val(i);
      gamma_vec[i] = tgamma(half_nu_dbl);
      digamma_vec[i] = digamma(half_nu_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_vec.val(n) == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return y_inv_dbl = 1.0 / y_dbl;
    const T_partials_return half_nu_dbl = 0.5 * nu_vec.val(n);
    const T_partials_return s_dbl = s_vec.val(n);
    const T_partials_return half_s2_overx_dbl = 0.5 * s_dbl * s_dbl * y_inv_dbl;
    const T_partials_return half_nu_s2_overx_dbl
        = 2.0 * half_nu_dbl * half_s2_overx_dbl;

    const T_partials_return Pn = gamma_p(half_nu_dbl, half_nu_s2_overx_dbl);
    const T_partials_return gamma_p_deriv
        = exp(-half_nu_s2_overx_dbl)
          * pow(half_nu_s2_overx_dbl, half_nu_dbl - 1) / tgamma(half_nu_dbl);

    P += log(Pn);

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[n]
          -= half_nu_s2_overx_dbl * y_inv_dbl * gamma_p_deriv / Pn;
    }
    if (!is_constant_all<T_dof>::value) {
      partials<1>(ops_partials)[n]
          -= (0.5
                  * grad_reg_inc_gamma(half_nu_dbl, half_nu_s2_overx_dbl,
                                       gamma_vec[n], digamma_vec[n])
              - half_s2_overx_dbl * gamma_p_deriv)
             / Pn;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<2>(ops_partials)[n]
          += 2.0 * half_nu_dbl * s_dbl * y_inv_dbl * gamma_p_deriv / Pn;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
