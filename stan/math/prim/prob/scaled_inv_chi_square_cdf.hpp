#ifndef STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_CDF_HPP
#define STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The CDF of a scaled inverse chi-squared density for y with the
 * specified degrees of freedom parameter and scale parameter.
 *
 * @param y A scalar variable.
 * @param nu Degrees of freedom.
 * @param s Scale parameter.
 * @throw std::domain_error if nu is not greater than 0
 * @throw std::domain_error if s is not greater than 0.
 * @throw std::domain_error if y is not greater than 0.
 * @tparam T_y Type of scalar.
 * @tparam T_dof Type of degrees of freedom.
 */
template <typename T_y, typename T_dof, typename T_scale>
return_type_t<T_y, T_dof, T_scale> scaled_inv_chi_square_cdf(const T_y& y,
                                                             const T_dof& nu,
                                                             const T_scale& s) {
  using T_partials_return = partials_return_t<T_y, T_dof, T_scale>;

  if (size_zero(y, nu, s)) {
    return 1.0;
  }

  static const char* function = "scaled_inv_chi_square_cdf";

  T_partials_return P(1.0);

  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_positive_finite(function, "Scale parameter", s);
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu, "Scale parameter",
                         s);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_dof> nu_vec(nu);
  scalar_seq_view<T_scale> s_vec(s);
  size_t N = max_size(y, nu, s);

  operands_and_partials<T_y, T_dof, T_scale> ops_partials(y, nu, s);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(y); i++) {
    if (value_of(y_vec[i]) == 0) {
      return ops_partials.build(0.0);
    }
  }

  using std::exp;
  using std::pow;

  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      gamma_vec(size(nu));
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_vec(size(nu));

  if (!is_constant_all<T_dof>::value) {
    for (size_t i = 0; i < size(nu); i++) {
      const T_partials_return half_nu_dbl = 0.5 * value_of(nu_vec[i]);
      gamma_vec[i] = tgamma(half_nu_dbl);
      digamma_vec[i] = digamma(half_nu_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) == INFTY) {
      continue;
    }

    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return y_inv_dbl = 1.0 / y_dbl;
    const T_partials_return half_nu_dbl = 0.5 * value_of(nu_vec[n]);
    const T_partials_return s_dbl = value_of(s_vec[n]);
    const T_partials_return half_s2_overx_dbl = 0.5 * s_dbl * s_dbl * y_inv_dbl;
    const T_partials_return half_nu_s2_overx_dbl
        = 2.0 * half_nu_dbl * half_s2_overx_dbl;

    const T_partials_return Pn = gamma_q(half_nu_dbl, half_nu_s2_overx_dbl);
    const T_partials_return gamma_p_deriv
        = exp(-half_nu_s2_overx_dbl)
          * pow(half_nu_s2_overx_dbl, half_nu_dbl - 1) / tgamma(half_nu_dbl);

    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += half_nu_s2_overx_dbl * y_inv_dbl * gamma_p_deriv / Pn;
    }

    if (!is_constant_all<T_dof>::value) {
      ops_partials.edge2_.partials_[n]
          += (0.5
                  * grad_reg_inc_gamma(half_nu_dbl, half_nu_s2_overx_dbl,
                                       gamma_vec[n], digamma_vec[n])
              - half_s2_overx_dbl * gamma_p_deriv)
             / Pn;
    }

    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += -2.0 * half_nu_dbl * s_dbl * y_inv_dbl * gamma_p_deriv / Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < size(y); ++n) {
      ops_partials.edge1_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_dof>::value) {
    for (size_t n = 0; n < size(nu); ++n) {
      ops_partials.edge2_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    for (size_t n = 0; n < size(s); ++n) {
      ops_partials.edge3_.partials_[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
