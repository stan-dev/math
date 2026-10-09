#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/prob/exp_mod_normal_lcdf.hpp>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_cdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_inv_scale&& lambda) {
  return exp(internal::exp_mod_normal_lcdf_impl<false>(
      "exp_mod_normal_cdf", std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(sigma), std::forward<T_inv_scale>(lambda)));
}

}  // namespace math
}  // namespace stan
#endif
