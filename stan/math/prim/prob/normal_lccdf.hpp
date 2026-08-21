#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP

#include <stan/math/prim/prob/normal_lcdf.hpp>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lccdf(const T_y& y,
                                                       const T_loc& mu,
                                                       const T_scale& sigma) {
  return normal_lcdf(-as_array_or_scalar(y), -as_array_or_scalar(mu), sigma);
}

}  // namespace math
}  // namespace stan
#endif
