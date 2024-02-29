#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_SQUARE_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_SQUARE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* inv_square_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_SQUARE\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_SQUARE\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Calculates `1 / (x*x)`
           *
           * <p><code>inv_square(x) = 1 / (x*x)</code>
           *
           * @param[in] x input
           * @return inverse square of the input
           *
           */
          double inv_square(double x) { return 1.0 / (x * x); }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
