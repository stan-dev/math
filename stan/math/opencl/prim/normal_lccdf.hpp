#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_LCCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/normal_lcdf.hpp>

namespace stan {
namespace math {

template <
    typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl> normal_lccdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  return internal::normal_lcdf_opencl_impl<true>("normal_lccdf(OpenCL)", y, mu,
                                                 sigma);
}

}  // namespace math
}  // namespace stan
#endif
#endif
