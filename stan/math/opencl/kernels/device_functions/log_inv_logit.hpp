#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* log_inv_logit_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_INV_LOGIT\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           *
           * Return the natural logarithm of the inverse logit of the
           * specified argument.
           *
           * @param x argument
           * @return the natural logarithm of the inverse logit of the
           *  specified argument
           */
          double log_inv_logit(double x) {
            if (x < 0.0) {
              return x - log1p_exp(x);  // prevent underflow
            }
            return -log1p_exp(-x);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
