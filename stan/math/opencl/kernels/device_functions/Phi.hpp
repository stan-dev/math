#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_PHI_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* phi_device_function
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
            // erfcx_cody_middle and erfcx_cody_tail come from the erfcx
            // device function, which the registration includes first.
            const double u = M_SQRT1_2 * x;
            const double a = fabs(u);
            if (a < 0.46875) {
              const double s = u * u;
              const double num
                  = (((1.85777706184603153e-1 * s + 3.16112374387056560) * s
                      + 1.13864154151050156e2)
                         * s
                     + 3.77485237685302021e2)
                        * s
                    + 3.20937758913846947e3;
              const double den
                  = (((s + 2.36012909523441209e1) * s + 2.44024637934444173e2)
                         * s
                     + 1.28261652607737228e3)
                        * s
                    + 2.84423683343917062e3;
              return 0.5 + 0.5 * (u * num / den);
            }
            const double tail
                = 0.5 * exp(-(0.5 * x) * x)
                  * (a >= 4.0 ? erfcx_cody_tail(a) : erfcx_cody_middle(a));
            return x < 0.0 ? tail : 1.0 - tail;
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
