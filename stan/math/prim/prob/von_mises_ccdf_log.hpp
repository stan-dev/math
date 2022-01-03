#ifndef STAN_MATH_PRIM_PROB_VON_MISES_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_CCDF_LOG_HPP

#include <stan/math/prim/prob/von_mises_lccdf.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>von_mises_lccdf</code>
 */
template <typename T_x, typename T_mu, typename T_k>
inline return_type_t<T_x, T_mu, T_k> von_mises_ccdf_log(const T_x& x,
                                                        const T_mu& mu,
                                                        const T_k& k) {
  return von_mises_lccdf<T_x, T_mu, T_k>(x, mu, k);
}

}  // namespace math
}  // namespace stan

#endif
