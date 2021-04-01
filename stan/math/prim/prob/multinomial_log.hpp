#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/multinomial_lpmf.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>multinomial_lpmf</code>
 */
template <bool propto, typename T_prob>
return_type_t<T_prob> multinomial_log(const std::vector<int>& ns,
                                      const T_prob& theta) {
  return multinomial_lpmf<propto>(ns, theta);
}

/** \ingroup multivar_dists
 * @deprecated use <code>multinomial_lpmf</code>
 */
template <typename T_prob>
return_type_t<T_prob> multinomial_log(const std::vector<int>& ns,
                                      const T_prob& theta) {
  return multinomial_lpmf<false>(ns, theta);
}

}  // namespace math
}  // namespace stan
#endif
