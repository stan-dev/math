#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_UTILS_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_UTILS_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/prob/normal_tail_utils.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

template <typename T>
struct exp_mod_normal_log_phi_terms {
  T log_term;
  T mills_ratio;
  T mills_excess;
};

template <typename T>
inline T exp_mod_normal_mills_excess(const T& x) {
  const T inv_x = 1.0 / x;
  const T inv_x_sq = inv_x * inv_x;
  T series = 706.0 - 8162.0 * inv_x_sq;
  series = -74.0 + inv_x_sq * series;
  series = 10.0 + inv_x_sq * series;
  series = -2.0 + inv_x_sq * series;
  return inv_x * (1.0 + inv_x_sq * series);
}

template <typename T>
inline exp_mod_normal_log_phi_terms<T> exp_mod_normal_exp_log_cdf_terms(
    const T& z, const T& a) {
  using std::log;
  const T erfc_arg = (a - z) * INV_SQRT_TWO;
  if (value_of_rec(erfc_arg) >= 5.0) {
    const T erfcx = erfcx_positive(erfc_arg);
    const T mills_excess = value_of_rec(erfc_arg) >= 20.0
                               ? exp_mod_normal_mills_excess(a - z)
                               : SQRT_TWO_OVER_SQRT_PI / erfcx - (a - z);
    return {-0.5 * z * z + LOG_HALF + log(erfcx), a - z + mills_excess,
            mills_excess};
  }
  const auto normal_terms = std_normal_lcdf_and_mills(z - a);
  return {0.5 * a * a - a * z + normal_terms.log_cdf, normal_terms.mills_ratio,
          normal_terms.mills_ratio - (a - z)};
}

template <typename T>
struct exp_mod_normal_log_cdf_terms {
  T log_cdf;
  T log_ccdf;
  T dz_log_cdf;
  T da_log_cdf;
  T a_da_log_cdf;
  T dz_log_ccdf;
  T da_log_ccdf;
};

template <typename T>
inline exp_mod_normal_log_cdf_terms<T> exp_mod_normal_cdf_terms(const T& z,
                                                                const T& a) {
  const auto z_terms = std_normal_lcdf_and_mills(z);
  const auto exp_cdf_terms = exp_mod_normal_exp_log_cdf_terms(z, a);
  using std::exp;
  using std::log;
  using std::log1p;
  if (value_of_rec(a) < 1e-8) {
    const T m0_factor = value_of_rec(z) < -4.0 ? exp_cdf_terms.mills_excess
                                               : z + z_terms.mills_ratio;
    const T m1_over_m0
        = 0.5 * (z * z + 1.0 + z * z_terms.mills_ratio) / m0_factor;
    if (value_of_rec(m0_factor) > 0.0
        && std::abs(value_of_rec(a * m1_over_m0)) < 1e-8) {
      const T remainder = 1.0 - a * m1_over_m0;
      const T log_m0 = z_terms.log_cdf + log(m0_factor);
      const T log_cdf = log(a) + log_m0 + log1p(-a * m1_over_m0);
      const T dm1_over_m0 = 1.0 - m1_over_m0 / m0_factor;
      const T dz_log_cdf = 1.0 / m0_factor - a * dm1_over_m0 / remainder;
      const T a_da_log_cdf = 1.0 - a * m1_over_m0 / remainder;
      const T da_log_cdf = a_da_log_cdf / a;
      const T log_ccdf = log1m_exp(log_cdf);
      const T cdf_to_ccdf = exp(log_cdf - log_ccdf);
      const T da_log_ccdf
          = -exp(log_m0 - log_ccdf) * (1.0 - 2.0 * a * m1_over_m0);
      return {log_cdf,    log_ccdf,     dz_log_cdf,
              da_log_cdf, a_da_log_cdf, -cdf_to_ccdf * dz_log_cdf,
              da_log_ccdf};
    }
  }

  const auto neg_z_terms = std_normal_lcdf_and_mills(-z);
  const T log_exp_cdf = exp_cdf_terms.log_term;
  const T log_cdf = log_diff_exp(z_terms.log_cdf, log_exp_cdf);
  const T log_ccdf = log_sum_exp(neg_z_terms.log_cdf, log_exp_cdf);

  const T cdf_weight = exp(z_terms.log_cdf - log_cdf);
  const T exp_cdf_weight = exp(log_exp_cdf - log_cdf);
  const T dz_log_cdf = cdf_weight * z_terms.mills_ratio
                       + exp_cdf_weight * (z - exp_cdf_terms.mills_excess);
  const T da_log_cdf = exp_cdf_weight * exp_cdf_terms.mills_excess;

  const T ccdf_weight = exp(neg_z_terms.log_cdf - log_ccdf);
  const T exp_ccdf_weight = exp(log_exp_cdf - log_ccdf);
  const T dz_log_ccdf = -ccdf_weight * neg_z_terms.mills_ratio
                        + exp_ccdf_weight * (-z + exp_cdf_terms.mills_excess);
  const T da_log_ccdf = -exp_ccdf_weight * exp_cdf_terms.mills_excess;

  return {log_cdf,        log_ccdf,    dz_log_cdf, da_log_cdf,
          a * da_log_cdf, dz_log_ccdf, da_log_ccdf};
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
