#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCCDF_HPP

#include <stan/math/prim/prob/std_normal_lcdf.hpp>

namespace stan {
namespace math {

template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_lccdf(T_y&& y) {
  return internal::std_normal_lcdf_impl<true>("std_normal_lccdf",
                                              std::forward<T_y>(y));
}

}  // namespace math
}  // namespace stan
#endif
