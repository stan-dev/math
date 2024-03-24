#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_PHI_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_PHI_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* inv_phi_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_PHI\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_PHI\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Return the inv_Phi function applied to the specified
           * argument.
           *
           * @return inv_Phi(x)
           */
          inline double inv_Phi(double p) {
            if (p < 8e-311) {
              return -INFINITY;
            }
            if (p == 1) {
              return INFINITY;
            }

            const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                 -2.759285104469687e+02, 1.383577518672690e+02,
                                 -3.066479806614716e+01, 2.506628277459239e+00};
            const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                 -1.556989798598866e+02, 6.680131188771972e+01,
                                 -1.328068155288572e+01};
            const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                 -2.400758277161838e+00, -2.549732539343734e+00,
                                 4.374664141464968e+00,  2.938163982698783e+00};
            const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                                 2.445134137142996e+00, 3.754408661907416e+00};

            const double p_low = 0.02425;
            const double p_high = 0.97575;

            double x;
            if ((p_low <= p) && (p <= p_high)) {
              double q = p - 0.5;
              double r = q * q;
              x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r
                   + a[5])
                  * q
                  / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4])
                         * r
                     + 1.0);
            } else if (p < p_low) {
              double q = sqrt(-2.0 * log(p));
              x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q
                   + c[5])
                  / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
            } else {
              double q = sqrt(-2.0 * log1m(p));
              x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q
                    + c[5])
                  / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
            }

            if (x < 37.6) {  // gradient blows up past here
              double e = Phi(x) - p;
              double u = e * sqrt(2.0 * M_PI) * exp(0.5 * x * x);
              x -= u / (1.0 + 0.5 * x * u);
            }

            return x;
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
