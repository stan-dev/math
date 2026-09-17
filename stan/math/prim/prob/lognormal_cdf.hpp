#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/prob/lognormal_lcdf.hpp>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> lognormal_cdf(T_y&& y, T_loc&& mu,
                                                        T_scale&& sigma) {
  return exp(internal::lognormal_lcdf_impl<false>(
      "lognormal_cdf", std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(sigma)));
}

}  // namespace math
}  // namespace stan
#endif
