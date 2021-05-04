#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_CCDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/bernoulli_lccdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>bernoulli_lccdf</code>
 */
template <typename T_n, typename T_prob>
return_type_t<T_prob> bernoulli_ccdf_log(const T_n& n, const T_prob& theta) {
  return bernoulli_lccdf<T_n, T_prob>(n, theta);
}
}  // namespace math
}  // namespace stan
#endif
