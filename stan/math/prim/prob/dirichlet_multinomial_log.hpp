#ifndef STAN_MATH_PRIM_PROB_DIRICHLET_MULTINOMIAL_LOG_HPP
#define STAN_MATH_PRIM_PROB_DIRICHLET_MULTINOMIAL_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/prob/dirichlet_multinomial_lpmf.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>dirichlet_multinomial_lpmf</code>
 */
template <bool propto, typename T_prior_size>
return_type_t<T_prior_size> dirichlet_multinomial_log(
    const std::vector<int>& ns, const T_prior_size& alpha) {
  return dirichlet_multinomial_lpmf<propto>(ns, alpha);
}

/** \ingroup multivar_dists
 * @deprecated use <code>dirichlet_multinomial_lpmf</code>
 */
template <typename T_prior_size>
return_type_t<T_prior_size> dirichlet_multinomial_log(
    const std::vector<int>& ns, const T_prior_size& alpha) {
  return dirichlet_multinomial_lpmf<false>(ns, alpha);
}

}  // namespace math
}  // namespace stan
#endif
