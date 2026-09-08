#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_LCCDF_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_LCCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/normal_lcdf.hpp>

namespace stan {
namespace math {
namespace internal {
constexpr char normal_lccdf_opencl_func[] = "normal_lccdf(OpenCL)";
}  // namespace internal

/** \ingroup opencl
 * Returns the normal log complementary cumulative distribution function
 * for the given location, and scale. If given containers of matching sizes
 * returns the log sum of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_loc_cl type of location
 * @tparam T_scale_cl type of scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <
    typename T_y_cl, typename T_loc_cl, typename T_scale_cl,
    require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                T_scale_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_y_cl, T_loc_cl, T_scale_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_scale_cl> normal_lccdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_scale_cl& sigma) {
  return normal_lcdf<internal::normal_lccdf_opencl_func>(-y, -mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
#endif
