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
          // Cody (1969) rationals, sharing the CPU kernel's coefficients.
          inline double std_normal_erf_small(double x) {
            const double a[] = {3.16112374387056560, 1.13864154151050156e2,
                                3.77485237685302021e2, 3.20937758913846947e3,
                                1.85777706184603153e-1};
            const double b[] = {2.36012909523441209e1, 2.44024637934444173e2,
                                1.28261652607737228e3, 2.84423683343917062e3};
            const double x2 = x * x;
            double numerator = a[4] * x2;
            double denominator = x2;
            for (int i = 0; i < 3; ++i) {
              numerator = (numerator + a[i]) * x2;
              denominator = (denominator + b[i]) * x2;
            }
            return x * (numerator + a[3]) / (denominator + b[3]);
          }

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

          /** erfcx(x) = exp(x^2) erfc(x) for x >= 0.46875. */
          inline double std_normal_erfcx(double x) {
            if (x > 4.0) {
              const double r = 1.0 / (x * x);
              return (1.0 + r * std_normal_lcdf_tail_correction(r))
                     * (0.5 * M_2_SQRTPI) / x;
            }
            const double c[] = {5.64188496988670089e-1, 8.88314979438837594,
                                6.61191906371416295e1,  2.98635138197400131e2,
                                8.81952221241769090e2,  1.71204761263407058e3,
                                2.05107837782607147e3,  1.23033935479799725e3,
                                2.15311535474403846e-8};
            const double d[] = {1.57449261107098347e1, 1.17693950891312499e2,
                                5.37181101862009858e2, 1.62138957456669019e3,
                                3.29079923573345963e3, 4.36261909014324716e3,
                                3.43936767414372164e3, 1.23033935480374942e3};
            double numerator = c[8] * x;
            double denominator = x;
            for (int i = 0; i < 7; ++i) {
              numerator = (numerator + c[i]) * x;
              denominator = (denominator + d[i]) * x;
            }
            return (numerator + c[7]) / (denominator + d[7]);
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
              const double e = std_normal_erf_small(s);
              return log1p(x < 0.0 ? -e : e) - M_LN2;
            }
            const double erfcx = std_normal_erfcx(s);
            if (x < 0.0) {
              return -(0.5 * x) * x - M_LN2 + log(erfcx);
            }
            return log1p(-0.5 * exp(-(0.5 * x) * x) * erfcx);
          }

          /** Slope phi(x) / Phi(x) in the original units. */
          inline double std_normal_lcdf_derivative(double x) {
            if (x <= -4.0 * M_SQRT2) {
              const double r = 2.0 / (x * x);
              return -x / (1.0 + r * std_normal_lcdf_tail_correction(r));
            }
            const double s = fabs(x) * M_SQRT1_2;
            if (s < 0.46875) {
              const double e = std_normal_erf_small(s);
              return (0.5 * M_SQRT2 * M_2_SQRTPI) * exp(-s * s)
                     / (x < 0.0 ? 1.0 - e : 1.0 + e);
            }
            const double erfcx = std_normal_erfcx(s);
            if (x < 0.0) {
              return (0.5 * M_SQRT2 * M_2_SQRTPI) / erfcx;
            }
            const double density = exp(-(0.5 * x) * x);
            return (0.5 * M_SQRT1_2 * M_2_SQRTPI) * density
                   / (1.0 - 0.5 * density * erfcx);
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
