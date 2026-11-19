#ifndef STAN_MATH_OPENCL_PRIM_LOGNORMAL_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_LOGNORMAL_CDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/opencl/prim/lognormal_lcdf.hpp>

namespace stan {
namespace math {

template <
    typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl> lognormal_cdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  return exp(internal::lognormal_lcdf_opencl_impl<false>(
      "lognormal_cdf(OpenCL)", y, mu, sigma));
}

}  // namespace math
}  // namespace stan
#endif
#endif
