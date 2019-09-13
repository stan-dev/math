#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_CCDF_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/prob/bernoulli_lccdf.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * @deprecated use <code>bernoulli_lccdf</code>
 */
template <typename T_n, typename T_prob>
inline auto bernoulli_ccdf_log(T_n&& n, T_prob&& theta) {
  return bernoulli_lccdf(std::forward<T_n>(n), std::forward<T_prob>(theta));
}
}  // namespace math
}  // namespace stan
#endif
