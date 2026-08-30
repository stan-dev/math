#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_LCCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/std_normal_lcdf.hpp>

namespace stan {
namespace math {
namespace internal {
constexpr char std_normal_lccdf_opencl_func[] = "std_normal_lccdf(OpenCL)";
}  // namespace internal

/** \ingroup opencl
 * Returns the log standard normal complementary cumulative distribution
 * function.
 *
 * @tparam T_y_cl type of scalar outcome
 * @param y (Sequence of) scalar(s).
 * @return The log of the product of densities.
 */
template <typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
inline return_type_t<T_y_cl> std_normal_lccdf(const T_y_cl& y) {
  return std_normal_lcdf<internal::std_normal_lccdf_opencl_func>(-y);
}

}  // namespace math
}  // namespace stan
#endif
#endif
