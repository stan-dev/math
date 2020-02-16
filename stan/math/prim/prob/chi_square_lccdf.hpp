#ifndef STAN_MATH_PRIM_PROB_CHI_SQUARE_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_CHI_SQUARE_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the chi square log complementary cumulative distribution
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
return_type_t<T_y, T_dof> chi_square_lccdf(const T_y& y, const T_dof& nu) {
  static const char* function = "chi_square_lccdf";
  using T_partials_return = partials_return_t<T_y, T_dof>;

  T_partials_return ccdf_log(0.0);

  if (size_zero(y, nu)) {
    return ccdf_log;
  }

  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_dof> nu_vec(nu);
  size_t N = max_size(y, nu);

  operands_and_partials<T_y, T_dof> ops_partials(y, nu);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(y); i++) {
    if (value_of(y_vec[i]) == 0) {
      return ops_partials.build(0.0);
    }
  }

  using std::exp;
  using std::log;
  using std::pow;

  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      gamma_vec(size(nu));
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_vec(size(nu));

  if (!is_constant_all<T_dof>::value) {
    for (size_t i = 0; i < size(nu); i++) {
      const T_partials_return alpha_dbl = value_of(nu_vec[i]) * 0.5;
      gamma_vec[i] = tgamma(alpha_dbl);
      digamma_vec[i] = digamma(alpha_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return alpha_dbl = value_of(nu_vec[n]) * 0.5;
    const T_partials_return beta_dbl = 0.5;

    const T_partials_return Pn = gamma_q(alpha_dbl, beta_dbl * y_dbl);

    ccdf_log += log(Pn);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= beta_dbl * exp(-beta_dbl * y_dbl)
                                          * pow(beta_dbl * y_dbl, alpha_dbl - 1)
                                          / tgamma(alpha_dbl) / Pn;
    }
    if (!is_constant_all<T_dof>::value) {
      ops_partials.edge2_.partials_[n]
          += 0.5
             * grad_reg_inc_gamma(alpha_dbl, beta_dbl * y_dbl, gamma_vec[n],
                                  digamma_vec[n])
             / Pn;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
