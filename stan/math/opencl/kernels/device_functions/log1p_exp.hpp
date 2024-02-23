#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1P_EXP_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1P_EXP_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* log1p_exp_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1P_EXP\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1P_EXP\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Calculates the log of 1 plus the exponential of the specified
           * value without overflow.
           *
           * <p><code>log1p_exp(x) = log(1+exp(x))</code>
           *
           * @param[in] a Argument.
           * @return natural logarithm of one plus the exponential of the
           * argument.
           */
          double log1p_exp(double a) {
            // prevents underflow
            return (a > 0 ? a : 0) + log1p(exp(-fabs(a)));
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
