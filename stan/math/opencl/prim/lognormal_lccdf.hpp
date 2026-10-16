#ifndef STAN_MATH_OPENCL_PRIM_LOGNORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_LOGNORMAL_LCCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/prim/lognormal_lcdf.hpp>

namespace stan {
namespace math {

template <
    typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl> lognormal_lccdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  return internal::lognormal_lcdf_opencl_impl<true>("lognormal_lccdf(OpenCL)",
                                                    y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
#endif
