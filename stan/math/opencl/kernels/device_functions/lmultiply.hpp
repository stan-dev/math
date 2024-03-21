#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LMULTIPLY_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LMULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* lmultiply_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LMULTIPLY\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LMULTIPLY\n" STRINGIFY(
          // \endcond
          /**
           * Return the first argument times the natural log of the second
           * argument if either argument is non-zero and 0 if both arguments
           * are 0.
           *
           * @param a first argument
           * @param b second argument
           * @return first argument times log of second argument
           */
          double lmultiply(double a, double b) {
            if (b == 0.0 && a == 0.0) {
              return 0.0;
            }
            return a * log(b);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
