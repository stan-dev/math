#ifndef STAN_MATH_OPENCL_PRIM_STD_NORMAL_CDF_HPP
#define STAN_MATH_OPENCL_PRIM_STD_NORMAL_CDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/opencl/prim/std_normal_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the standard normal cumulative distribution function. If given a
 * container, returns the product of probabilities.
 *
 * @tparam T_y_cl type of scalar outcome
 * @param y (Sequence of) scalar(s).
 * @return The product of cumulative probabilities.
 */
template <typename T_y_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl>* = nullptr>
inline return_type_t<T_y_cl> std_normal_cdf(const T_y_cl& y) {
  return exp(internal::std_normal_lcdf_opencl_impl<false>(
      "std_normal_cdf(OpenCL)", y));
}

}  // namespace math
}  // namespace stan
#endif
#endif
