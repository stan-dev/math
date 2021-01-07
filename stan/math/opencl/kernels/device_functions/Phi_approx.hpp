#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_APPROX_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_APPROX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* phi_approx_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_APPROX\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_APPROX\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Return the Phi_approx function applied to the specified
           * argument.
           *
           * @return Phi_approx(x)
           */
          inline double Phi_approx(double x) {
            return inv_logit(0.07056 * pow(x, 3.0) + 1.5976 * x);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
