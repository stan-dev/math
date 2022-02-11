#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* lgamma_stirling_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING\n" STRINGIFY(
          // \endcond
          /**
           * Return the Stirling approximation to the lgamma function.
           *
             \f[
             \mbox{lgamma_stirling}(x) =
              \frac{1}{2} \log(2\pi) + (x-\frac{1}{2})*\log(x) - x
             \f]
           *
           * @param x value
           * @return Stirling's approximation to lgamma(x).
           */
          double lgamma_stirling(double x) {
            return 0.5 * log(2.0 * M_PI) + (x - 0.5) * log(x) - x;
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
