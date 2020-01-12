#ifndef STAN_MATH_PRIM_PROB_INV_CHI_SQUARE_LPDF_HPP
#define STAN_MATH_PRIM_PROB_INV_CHI_SQUARE_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of an inverse chi-squared density for y with the specified
 * degrees of freedom parameter.
 * The degrees of freedom prarameter must be greater than 0.
 * y must be greater than 0.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{Inv-}}\chi^2_\nu \\
 \log (p (y \, |\, \nu)) &=& \log \left( \frac{2^{-\nu / 2}}{\Gamma (\nu / 2)}
 y^{- (\nu / 2 + 1)} \exp^{-1 / (2y)} \right) \\
 &=& - \frac{\nu}{2} \log(2) - \log (\Gamma (\nu / 2)) - (\frac{\nu}{2} + 1)
 \log(y) - \frac{1}{2y} \\ & & \mathrm{ where } \; y > 0 \f}
 * @param y A scalar variable.
 * @param nu Degrees of freedom.
 * @throw std::domain_error if nu is not greater than or equal to 0
 * @throw std::domain_error if y is not greater than or equal to 0.
 * @tparam T_y Type of scalar.
 * @tparam T_dof Type of degrees of freedom.
 */
template <bool propto, typename T_y, typename T_dof>
return_type_t<T_y, T_dof> inv_chi_square_lpdf(const T_y& y, const T_dof& nu) {
  static const char* function = "inv_chi_square_lpdf";
  using T_partials_return = partials_return_t<T_y, T_dof>;

  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_not_nan(function, "Random variable", y);
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);
  if (size_zero(y, nu)) {
    return 0;
  }

  T_partials_return logp(0);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_dof> nu_vec(nu);
  size_t N = max_size(y, nu);

  for (size_t n = 0; n < size(y); n++) {
    if (value_of(y_vec[n]) <= 0) {
      return LOG_ZERO;
    }
  }

  using std::log;

  VectorBuilder<include_summand<propto, T_y, T_dof>::value, T_partials_return,
                T_y>
      log_y(size(y));
  for (size_t i = 0; i < size(y); i++) {
    if (include_summand<propto, T_y, T_dof>::value) {
      log_y[i] = log(value_of(y_vec[i]));
    }
  }

  VectorBuilder<include_summand<propto, T_y>::value, T_partials_return, T_y>
      inv_y(size(y));
  for (size_t i = 0; i < size(y); i++) {
    if (include_summand<propto, T_y>::value) {
      inv_y[i] = 1.0 / value_of(y_vec[i]);
    }
  }

  VectorBuilder<include_summand<propto, T_dof>::value, T_partials_return, T_dof>
      lgamma_half_nu(size(nu));
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_half_nu_over_two(size(nu));
  for (size_t i = 0; i < size(nu); i++) {
    T_partials_return half_nu = 0.5 * value_of(nu_vec[i]);
    if (include_summand<propto, T_dof>::value) {
      lgamma_half_nu[i] = lgamma(half_nu);
    }
    if (!is_constant_all<T_dof>::value) {
      digamma_half_nu_over_two[i] = digamma(half_nu) * 0.5;
    }
  }

  operands_and_partials<T_y, T_dof> ops_partials(y, nu);
  for (size_t n = 0; n < N; n++) {
    const T_partials_return nu_dbl = value_of(nu_vec[n]);
    const T_partials_return half_nu = 0.5 * nu_dbl;

    if (include_summand<propto, T_dof>::value) {
      logp -= nu_dbl * HALF_LOG_TWO + lgamma_half_nu[n];
    }
    if (include_summand<propto, T_y, T_dof>::value) {
      logp -= (half_nu + 1.0) * log_y[n];
    }
    if (include_summand<propto, T_y>::value) {
      logp -= 0.5 * inv_y[n];
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += -(half_nu + 1.0) * inv_y[n] + 0.5 * inv_y[n] * inv_y[n];
    }
    if (!is_constant_all<T_dof>::value) {
      ops_partials.edge2_.partials_[n]
          -= HALF_LOG_TWO + digamma_half_nu_over_two[n] + 0.5 * log_y[n];
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_dof>
inline return_type_t<T_y, T_dof> inv_chi_square_lpdf(const T_y& y,
                                                     const T_dof& nu) {
  return inv_chi_square_lpdf<false>(y, nu);
}

}  // namespace math
}  // namespace stan
#endif
