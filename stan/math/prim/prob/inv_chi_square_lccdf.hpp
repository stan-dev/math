#ifndef STAN_MATH_PRIM_PROB_INV_CHI_SQUARE_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_INV_CHI_SQUARE_LCCDF_HPP

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

/** \ingroup prob_dists
 * Returns the inverse chi square log complementary cumulative distribution
 * function for the given variate and degrees of freedom. If given
 * containers of matching sizes, returns the log sum of probabilities.
 *
 * @tparam T_y type of scalar parameter
 * @tparam T_dof type of degrees of freedom parameter
 * @param y scalar parameter
 * @param nu degrees of freedom parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is negative or nu is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_dof>
return_type_t<T_y, T_dof> inv_chi_square_lccdf(const T_y& y, const T_dof& nu) {
  using T_partials_return = partials_return_t<T_y, T_dof>;
  using std::exp;
  using std::log;
  using std::pow;
  static const char* function = "inv_chi_square_lccdf";
  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);

  if (size_zero(y, nu)) {
    return 0;
  }

  T_partials_return P(0.0);
  operands_and_partials<T_y, T_dof> ops_partials(y, nu);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_dof> nu_vec(nu);
  size_t size_y = stan::math::size(y);
  size_t size_nu = stan::math::size(nu);
  size_t N = max_size(y, nu);

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

  VectorBuilder<true, T_partials_return, T_y> half_y_inv(size_y);
  for (size_t i = 0; i < size_y; i++) {
    half_y_inv[i] = 0.5 * inv(value_of(y_vec[i]));
  }

  VectorBuilder<true, T_partials_return, T_dof> half_nu(size_nu);
  VectorBuilder<!is_constant_all<T_y, T_dof>::value, T_partials_return, T_dof>
      tgamma_vec(size_nu);
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_vec(size_nu);
  for (size_t i = 0; i < size_nu; i++) {
    const T_partials_return half_nu_dbl = 0.5 * value_of(nu_vec[i]);
    half_nu[i] = half_nu_dbl;
    if (!is_constant_all<T_y, T_dof>::value) {
      tgamma_vec[i] = tgamma(half_nu[i]);
    }
    if (!is_constant_all<T_dof>::value) {
      digamma_vec[i] = digamma(half_nu_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return Pn = gamma_p(half_nu[n], half_y_inv[n]);
    P += log(Pn);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          -= 2 * square(half_y_inv[n]) * exp(-half_y_inv[n])
             * pow(half_y_inv[n], half_nu[n] - 1) / (tgamma_vec[n] * Pn);
    }
    if (!is_constant_all<T_dof>::value) {
      ops_partials.edge2_.partials_[n]
          -= 0.5
             * grad_reg_inc_gamma(half_nu[n], half_y_inv[n], tgamma_vec[n],
                                  digamma_vec[n])
             / Pn;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
