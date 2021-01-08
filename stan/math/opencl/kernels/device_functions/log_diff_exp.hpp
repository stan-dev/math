#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_DIFF_EXP_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_DIFF_EXP_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* log_diff_exp_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_DIFF_EXP\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_DIFF_EXP\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           *
           * Returns the natural logarithm of the difference of the
           * inverse logits of the specified arguments.
           *
           * @param x argument
           * @return the natural logarithm of the inverse logit of the
           *  specified argument
           */
          double log_diff_exp(double x, double y) {
            if (x <= y) {
              return (x < INFINITY && x == y) ? -INFINITY : NAN;
            }
            return x + log1m_exp(y - x);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
