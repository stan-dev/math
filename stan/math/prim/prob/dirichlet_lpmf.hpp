#ifndef STAN_MATH_PRIM_PROB_DIRICHLET_LPMF_HPP
#define STAN_MATH_PRIM_PROB_DIRICHLET_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/dirichlet_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * @deprecated use <code>dirichlet_lpdf</code>
 */
template <bool propto, typename T_prob, typename T_prior_size>
return_type_t<T_prob, T_prior_size> dirichlet_lpmf(const T_prob& theta,
                                                   const T_prior_size& alpha) {
  return dirichlet_lpdf<propto, T_prob, T_prior_size>(theta, alpha);
}

template <typename T_prob, typename T_prior_size>
return_type_t<T_prob, T_prior_size> dirichlet_lpmf(const T_prob& theta,
                                                   const T_prior_size& alpha) {
  return dirichlet_lpdf<T_prob, T_prior_size>(theta, alpha);
}

}  // namespace math
}  // namespace stan
#endif
