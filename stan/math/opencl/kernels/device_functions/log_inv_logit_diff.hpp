#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT_DIFF_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT_DIFF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* log_inv_logit_diff_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT_DIFF\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT_"
      "DIFF\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           *
           * Returns the natural logarithm of the difference of the
           * inverse logits of the specified arguments.
           *
           * @param x first argument
           * @param y second argument
           * @return the natural logarithm of the difference of the
           * inverse logits of the specified argument.
           */
          double log_inv_logit_diff(double x, double y) {
            return x - log1p_exp(x) + log1m_exp(y - x) - log1p_exp(y);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
