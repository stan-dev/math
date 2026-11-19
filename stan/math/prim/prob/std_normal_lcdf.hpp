#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <utility>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * @brief Calculates the log of the cdf of the standard normal distribution
 *
 * Shares the scalar value and slope calculation with normal_lcdf through
 * std_normal_lcdf_impl.hpp.
 *
 * @tparam T_y A vector or scalar type for the random variable.
 * @param y (Sequence of) scalar(s).
 * @return The log of the standard normal cdf evaluated at the specified
 *   argument. If given a container, the log of the product of the cdfs.
 */
template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_lcdf(T_y&& y) {
  return internal::std_normal_lcdf_impl<false>("std_normal_lcdf",
                                               std::forward<T_y>(y));
}

}  // namespace math
}  // namespace stan
#endif
