#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP

#include <stan/math/prim/prob/normal_lcdf.hpp>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lccdf(T_y&& y, T_loc&& mu,
                                                       T_scale&& sigma) {
  return internal::normal_lcdf_impl<true>("normal_lccdf", std::forward<T_y>(y),
                                          std::forward<T_loc>(mu),
                                          std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
