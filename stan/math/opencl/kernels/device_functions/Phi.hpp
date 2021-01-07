#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* phi_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Return the Phi function applied to the specified
           * argument.
           *
           * @return Phi(x)
           */
          inline double Phi(double x) {
            if (x < -37.5) {
              return 0;
            } else if (x < -5.0) {
              return 0.5 * erfc(-1.0 / sqrt(2.0) * x);
            } else if (x > 8.25) {
              return 1;
            } else {
              return 0.5 * (1.0 + erf(1.0 / sqrt(2.0) * x));
            }
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
