#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/std_normal_lcdf.hpp>

namespace stan {
namespace math {

template <typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
inline return_type_t<T_y_cl> std_normal_lccdf(const T_y_cl& y) {
  return internal::std_normal_lcdf_opencl_impl<true>("std_normal_lccdf(OpenCL)",
                                                     y);
}

}  // namespace math
}  // namespace stan
#endif
#endif
