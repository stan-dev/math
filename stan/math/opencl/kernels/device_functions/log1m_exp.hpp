#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_EXP_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_EXP_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* log1m_exp_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_EXP\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_EXP\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Calculates the natural logarithm of one minus the exponential
           * of the specified value without overflow,
           *
           * <p><code>log1m_exp(x) = log(1-exp(x))</code>
           *
           * This function is only defined for x < 0
           *
           * @param[in] a Argument.
           * @return natural logarithm of one minus the exponential of the
           * argument.
           *
           */
          double log1m_exp(double a) {
            if (a > -0.693147)
              return log(-expm1(a));  // 0.693147 ~= log(2)
            else
              return log1p(-exp(a));
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
