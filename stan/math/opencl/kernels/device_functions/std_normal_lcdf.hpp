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
          // Cody tail correction, shared with the CPU kernel's coefficients.
          inline double std_normal_lcdf_tail_correction(double r) {
            const double p[]
                = {0.000658749161529837803157, 0.0160837851487422766278,
                   0.125781726111229246204,    0.360344899949804439429,
                   0.305326634961232344035,    0.0163153871373020978498};
            const double q[]
                = {-0.00233520497626869185443, -0.0605183413124413191178,
                   -0.527905102951428412248,   -1.87295284992346047209,
                   -2.56852019228982242072,    -1.0};
            double numerator = p[5] * r + p[4];
            double denominator = q[5] * r + q[4];
            for (int i = 3; i >= 0; --i) {
              numerator = numerator * r + p[i];
              denominator = denominator * r + q[i];
            }
            return (numerator / denominator) / (0.5 * M_2_SQRTPI);
          }

          /** Log Phi(x), with the original, unscaled argument. */
          inline double std_normal_lcdf_impl(double x) {
            if (x > 0.0) {
              return log1p(-0.5 * erfc(x * M_SQRT1_2));
            }
            if (x > -4.0 * M_SQRT2) {
              return log(0.5 * erfc(-x * M_SQRT1_2));
            }
            const double inv_a = -1.0 / x;
            const double r = 2 * inv_a * inv_a;
            return -(0.5 * x) * x - log(-x) - 0.91893853320467274178
                   + log1p(r * std_normal_lcdf_tail_correction(r));
          }

          /** Slope in the original units; scaling it up first can overflow. */
          inline double std_normal_lcdf_derivative(double x) {
            if (x <= -4.0 * M_SQRT2) {
              const double inv_a = -1.0 / x;
              const double r = 2 * inv_a * inv_a;
              const double correction = std_normal_lcdf_tail_correction(r);
              return -x - 2 * correction * inv_a / (1 + r * correction);
            }
            return (0.5 * M_SQRT2 * M_2_SQRTPI) * exp(-(0.5 * x) * x)
                   / erfc(-x * M_SQRT1_2);
          }

          // Compatibility for callers parameterized in units of sqrt(2).
          inline double std_normal_lcdf_scaled_impl(double x) {
            return std_normal_lcdf_impl(x * M_SQRT2);
          } inline double std_normal_lcdf_dscaled_impl(double x) {
            return M_SQRT2 * std_normal_lcdf_derivative(x * M_SQRT2);
          }) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
