#ifndef STAN_MATH_PRIM_PROB_BETA_BINOMIAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_PROB_BETA_BINOMIAL_CDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/beta_binomial_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @deprecated use <code>beta_binomial_lcdf</code>
 */
template <typename T_n, typename T_N, typename T_size1, typename T_size2>
return_type_t<T_size1, T_size2> beta_binomial_cdf_log(const T_n& n,
                                                      const T_N& N,
                                                      const T_size1& alpha,
                                                      const T_size2& beta) {
  return beta_binomial_lcdf<T_n, T_N, T_size1, T_size2>(n, N, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
