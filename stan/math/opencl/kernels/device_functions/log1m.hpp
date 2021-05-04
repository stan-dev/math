#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* log1m_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Calculates the natural logarithm of one minus
           * the specified value.
           *
           * <p><code>log1m(x) = log1p(-x)</code>
           *
           * @param[in] a Argument.
           * @return natural logarithm of one minus the argument.
           *
           */
          double log1m(double a) { return log1p(-a); }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
