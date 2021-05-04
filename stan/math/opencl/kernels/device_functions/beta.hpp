#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* beta_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Return the beta function applied to the specified
           * arguments.
           *
           * @param a First value
           * @param b Second value
           * @return beta(a, b)
           */
          double beta(double a, double b) {
            return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
