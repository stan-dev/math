#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_STD_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_STD_NORMAL_LCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* std_normal_lcdf_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_STD_NORMAL_LCDF\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_STD_NORMAL_LCDF\n" STRINGIFY(
          // The erfcx pieces come from the erfcx device function, which the
          // registration includes first.
          inline double std_normal_lcdf_tail_correction(double r) {
            return erfcx_tail_correction(r) / (0.5 * M_2_SQRTPI);
          }

          /** erfcx(x) = exp(x^2) erfc(x) for x >= 0.46875. */
          inline double std_normal_erfcx(double x) {
            return x >= 4.0 ? erfcx_cody_tail(x) : erfcx_cody_middle(x);
          }

          /** Log Phi(x), with the original, unscaled argument. */
          inline double std_normal_lcdf_impl(double x) {
            if (x <= -4.0 * M_SQRT2) {
              const double r = 2.0 / (x * x);
              return -(0.5 * x) * x - log(-x) - 0.91893853320467274178
                     + log1p(r * std_normal_lcdf_tail_correction(r));
            }
            const double s = fabs(x) * M_SQRT1_2;
            if (s < 0.46875) {
              // Phi(x) = exp(-x^2 / 2) erfcx(-x / sqrt(2)) / 2.
              return -M_LN2 - (0.5 * x) * x + log(erfcx_small(-x * M_SQRT1_2));
            }
            const double erfcx_s = std_normal_erfcx(s);
            if (x < 0.0) {
              return -(0.5 * x) * x - M_LN2 + log(erfcx_s);
            }
            return log1p(-0.5 * exp(-(0.5 * x) * x) * erfcx_s);
          }

          /** Slope phi(x) / Phi(x) in the original units. */
          inline double std_normal_lcdf_derivative(double x) {
            if (x <= -4.0 * M_SQRT2) {
              const double r = 2.0 / (x * x);
              return -x / (1.0 + r * std_normal_lcdf_tail_correction(r));
            }
            const double s = fabs(x) * M_SQRT1_2;
            if (s < 0.46875) {
              return (0.5 * M_SQRT2 * M_2_SQRTPI) / erfcx_small(-x * M_SQRT1_2);
            }
            const double erfcx_s = std_normal_erfcx(s);
            if (x < 0.0) {
              return (0.5 * M_SQRT2 * M_2_SQRTPI) / erfcx_s;
            }
            const double density = exp(-(0.5 * x) * x);
            return (0.5 * M_SQRT1_2 * M_2_SQRTPI) * density
                   / (1.0 - 0.5 * density * erfcx_s);
          }) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
