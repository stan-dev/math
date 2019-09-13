#ifndef STAN_MATH_PRIM_SCAL_PROB_CHI_SQUARE_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CHI_SQUARE_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/tgamma.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_gamma.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Returns the chi square log cumulative distribution function for the given
 * variate and degrees of freedom. If given containers of matching sizes,
 * returns the log sum of probabilities.
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
inline auto chi_square_lcdf(const T_y& y, const T_dof& nu) {
  using T_partials = partials_return_t<T_y, T_dof>;
  T_partials cdf_log(0.0);
  using T_return = return_type_t<T_y, T_dof>;
  using std::exp;
  using std::log;
  using std::pow;

  static const char* function = "chi_square_lcdf";
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_dof> nu_vec(nu);
  const size_t N = max_size(y, nu);

  operands_and_partials<T_y, T_dof> ops_partials(y, nu);
  if (size_zero(y, nu)) {
    return ops_partials.build(cdf_log);
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials y_dbl = value_of(y_vec[n]);
    const T_partials alpha_dbl = value_of(nu_vec[n]) * 0.5;
    const T_partials beta_dbl = 0.5;
    const T_partials Pn = gamma_p(alpha_dbl, beta_dbl * y_dbl);
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_dbl == std::numeric_limits<double>::infinity()) {
      return ops_partials.build(T_partials(0.0));
    }
    if (y_dbl == 0) {
      return ops_partials.build(negative_infinity());
    }

    cdf_log += log(Pn);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += beta_dbl * exp(-beta_dbl * y_dbl)
                                          * pow(beta_dbl * y_dbl, alpha_dbl - 1)
                                          / tgamma(alpha_dbl) / Pn;
    }
    if (!is_constant_all<T_dof>::value) {
      const T_partials gamma_dbl = tgamma(alpha_dbl);
      const T_partials digamma_dbl = digamma(alpha_dbl);
      ops_partials.edge2_.partials_[n]
          -= 0.5
             * grad_reg_inc_gamma(alpha_dbl, beta_dbl * y_dbl, gamma_dbl,
                                  digamma_dbl)
             / Pn;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
