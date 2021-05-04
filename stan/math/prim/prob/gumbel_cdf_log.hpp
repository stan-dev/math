#ifndef STAN_MATH_PRIM_PROB_GUMBEL_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_GUMBEL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/gumbel_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>gumbel_lcdf</code>
 */
template <typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> gumbel_cdf_log(const T_y& y, const T_loc& mu,
                                                  const T_scale& beta) {
  return gumbel_lcdf<T_y, T_loc, T_scale>(y, mu, beta);
}

}  // namespace math
}  // namespace stan
#endif
