#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/prob/std_normal_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the standard normal cumulative distribution function
 * for the given variate.
 *
 * \f$\Phi(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\inf}^x e^{-t^2/2} dt\f$.
 *
 * Evaluated as the exponential of the log cdf so the tails and gradients
 * share the log kernel instead of hard cutoffs.
 *
 * @tparam T_y type of y
 * @param y scalar variate
 * @return The standard normal cdf evaluated at the specified argument.
 */
template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_cdf(T_y&& y) {
  return exp(internal::std_normal_lcdf_impl<false>("std_normal_cdf",
                                                   std::forward<T_y>(y)));
}

}  // namespace math
}  // namespace stan
#endif
