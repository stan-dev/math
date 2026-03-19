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
          /** \ingroup opencl_kernels
           * Return the log standard normal cumulative distribution function
           * evaluated from the scaled input `x / sqrt(2)`.
           *
           * @param scaled_y input scaled by `1 / sqrt(2)`
           * @return log(Phi(x))
           */
          inline double std_normal_lcdf_scaled_impl(double scaled_y) {
            double lcdf_n;
            if (scaled_y > 0.0) {
              // CDF(x) = 1/2 + 1/2 erf(x) = 1 - 1/2 erfc(x)
              lcdf_n = log1p(-0.5 * erfc(scaled_y));
              if (isnan(lcdf_n)) {
                lcdf_n = 0;
              }
            } else if (scaled_y > -20.0) {
              // CDF(x) = 1/2 - 1/2 erf(-x) = 1/2 erfc(-x)
              lcdf_n = log(erfc(-scaled_y)) - M_LN2;
            } else if (10.0 * log(fabs(scaled_y)) < log(DBL_MAX)) {
              // Need direct approximation once erfc(-x) underflows.
              const double x2 = scaled_y * scaled_y;
              const double x4 = pow(scaled_y, 4);
              const double x6 = pow(scaled_y, 6);
              const double x8 = pow(scaled_y, 8);
              const double x10 = pow(scaled_y, 10);
              const double temp_p = 0.000658749161529837803157
                                    + 0.0160837851487422766278 / x2
                                    + 0.125781726111229246204 / x4
                                    + 0.360344899949804439429 / x6
                                    + 0.305326634961232344035 / x8
                                    + 0.0163153871373020978498 / x10;
              const double temp_q
                  = -0.00233520497626869185443 - 0.0605183413124413191178 / x2
                    - 0.527905102951428412248 / x4 - 1.87295284992346047209 / x6
                    - 2.56852019228982242072 / x8 - 1.0 / x10;
              lcdf_n = log(0.5 * M_2_SQRTPI + (temp_p / temp_q) / x2) - M_LN2
                       - log(-scaled_y) - x2;
            } else {
              lcdf_n = -INFINITY;
            }
            return lcdf_n;
          }

          /** \ingroup opencl_kernels
           * Return the derivative of log standard normal cumulative
           * distribution function with respect to the scaled input
           * `x / sqrt(2)`.
           *
           * @param scaled_y input scaled by `1 / sqrt(2)`
           * @return d / d(scaled_y) log(Phi(x))
           */
          inline double std_normal_lcdf_dscaled_impl(double scaled_y) {
            double dnlcdf = 0.0;
            double t = 0.0;
            double t2 = 0.0;
            double t4 = 0.0;
            const double x2 = scaled_y * scaled_y;

            if (scaled_y > 2.9) {
              t = 1.0 / (1.0 + 0.3275911 * scaled_y);
              t2 = t * t;
              t4 = pow(t, 4);
              dnlcdf = 0.5 * M_2_SQRTPI
                       / (exp(x2) - 0.254829592 + 0.284496736 * t
                          - 1.421413741 * t2 + 1.453152027 * t2 * t
                          - 1.061405429 * t4);
            } else if (scaled_y > 2.5) {
              t = scaled_y - 2.7;
              t2 = t * t;
              t4 = pow(t, 4);
              dnlcdf = 0.0003849882382 - 0.002079084702 * t
                       + 0.005229340880 * t2 - 0.008029540137 * t2 * t
                       + 0.008232190507 * t4 - 0.005692364250 * t4 * t
                       + 0.002399496363 * pow(t, 6);
            } else if (scaled_y > 2.1) {
              t = scaled_y - 2.3;
              t2 = t * t;
              t4 = pow(t, 4);
              dnlcdf = 0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
                       - 0.03326906904 * t2 * t + 0.02482478940 * t4
                       - 0.009883071924 * t4 * t - 0.0002771362254 * pow(t, 6);
            } else if (scaled_y > 1.5) {
              t = scaled_y - 1.85;
              t2 = t * t;
              t4 = pow(t, 4);
              dnlcdf = 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
                       - 0.09274533184 * t2 * t + 0.03543327418 * t4
                       + 0.005644855518 * t4 * t - 0.01111434424 * pow(t, 6);
            } else if (scaled_y > 0.8) {
              t = scaled_y - 1.15;
              t2 = t * t;
              t4 = pow(t, 4);
              dnlcdf = 0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
                       - 0.09748053605 * t2 * t - 0.04347986191 * t4
                       + 0.02182506378 * t4 * t + 0.01074751427 * pow(t, 6);
            } else if (scaled_y > 0.1) {
              t = scaled_y - 0.45;
              t2 = t * t;
              t4 = pow(t, 4);
              dnlcdf = 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
                       + 0.04700850676 * t2 * t - 0.03478651979 * t4
                       - 0.01772675404 * t4 * t + 0.0006577254811 * pow(t, 6);
            } else if (10.0 * log(fabs(scaled_y)) < log(DBL_MAX)) {
              t = 1.0 / (1.0 - 0.3275911 * scaled_y);
              t2 = t * t;
              t4 = pow(t, 4);
              dnlcdf
                  = M_2_SQRTPI
                    / (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
                       - 1.453152027 * t4 + 1.061405429 * t4 * t);
              if (scaled_y < -29.0) {
                dnlcdf += 0.0015065154280332 * x2
                          - 0.3993154819705530 * scaled_y - 4.2919418242931700;
              } else if (scaled_y < -17.0) {
                dnlcdf += 0.0001263257217272 * x2 * scaled_y
                          + 0.0123586859488623 * x2
                          - 0.0860505264736028 * scaled_y - 1.252783383752970;
              } else if (scaled_y < -7.0) {
                dnlcdf += 0.000471585349920831 * x2 * scaled_y
                          + 0.0296839305424034 * x2
                          + 0.207402143352332 * scaled_y + 0.425316974683324;
              } else if (scaled_y < -3.9) {
                dnlcdf += -0.0006972280656443 * x2 * scaled_y
                          + 0.0068218494628567 * x2
                          + 0.0585761964460277 * scaled_y + 0.1034397670201370;
              } else if (scaled_y < -2.1) {
                dnlcdf += -0.0018742199480885 * x2 * scaled_y
                          - 0.0097119598291202 * x2
                          - 0.0170137970924080 * scaled_y - 0.0100428567412041;
              }
            } else {
              dnlcdf = INFINITY;
            }

            return dnlcdf;
          }) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
