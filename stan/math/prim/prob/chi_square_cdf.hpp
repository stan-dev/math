#ifndef STAN_MATH_PRIM_PROB_CHI_SQUARE_CDF_HPP
#define STAN_MATH_PRIM_PROB_CHI_SQUARE_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the chi square cumulative distribution function for the given
 * variate and degrees of freedom. If given containers of matching sizes,
 * returns the product of probabilities.
 *
 * @tparam T_y type of scalar parameter
 * @tparam T_dof type of degrees of freedom parameter
 * @param y scalar parameter
 * @param nu degrees of freedom parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if y is negative or nu is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_dof>
return_type_t<T_y, T_dof> chi_square_cdf(const T_y& y, const T_dof& nu) {
  static const char* function = "chi_square_cdf";
  using T_partials_return = partials_return_t<T_y, T_dof>;

  T_partials_return cdf(1.0);

  if (size_zero(y, nu)) {
    return cdf;
  }

  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_dof> nu_vec(nu);
  size_t size_nu = stan::math::size(nu);
  size_t N = max_size(y, nu);

  operands_and_partials<T_y, T_dof> ops_partials(y, nu);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (value_of(y_vec[i]) == 0) {
      return ops_partials.build(0.0);
    }
  }

  using std::exp;
  using std::pow;

  VectorBuilder<!is_constant_all<T_y, T_dof>::value, T_partials_return, T_dof>
      tgamma_vec(size_nu);
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_vec(size_nu);
  if (!is_constant_all<T_y, T_dof>::value) {
    for (size_t i = 0; i < size_nu; i++) {
      const T_partials_return half_nu_dbl = 0.5 * value_of(nu_vec[i]);
      tgamma_vec[i] = tgamma(half_nu_dbl);
      if (!is_constant_all<T_dof>::value) {
        digamma_vec[i] = digamma(half_nu_dbl);
      }
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) == INFTY) {
      continue;
    }

    const T_partials_return half_nu_dbl = 0.5 * value_of(nu_vec[n]);
    const T_partials_return half_y_dbl = 0.5 * value_of(y_vec[n]);
    const T_partials_return Pn = gamma_p(half_nu_dbl, half_y_dbl);

    cdf *= Pn;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += 0.5 * exp(-half_y_dbl)
                                          * pow(half_y_dbl, half_nu_dbl - 1)
                                          / (tgamma_vec[n] * Pn);
    }
    if (!is_constant_all<T_dof>::value) {
      ops_partials.edge2_.partials_[n]
          -= 0.5
             * grad_reg_inc_gamma(half_nu_dbl, half_y_dbl, tgamma_vec[n],
                                  digamma_vec[n])
             / Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      ops_partials.edge1_.partials_[n] *= cdf;
    }
  }
  if (!is_constant_all<T_dof>::value) {
    for (size_t n = 0; n < size_nu; ++n) {
      ops_partials.edge2_.partials_[n] *= cdf;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
