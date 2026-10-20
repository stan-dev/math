#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_EXP_MOD_NORMAL_LCDF_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_EXP_MOD_NORMAL_LCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* exp_mod_normal_lcdf_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_EXP_MOD_NORMAL_LCDF\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_EXP_MOD_NORMAL_"
      "LCDF\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * For x >= 4, erfcx(x) = (1 + h(x)) / (sqrt(pi) x).
           *
           * Mirrors stan/math/prim/prob/exp_mod_normal_lcdf.hpp, which
           * documents the formulas; the two must be kept in step. The
           * erfcx functions come from the erfcx device function, which the
           * registration includes first.
           *
           * @param x argument, x >= 4
           * @return h(x)
           */
          double exp_mod_normal_erfcx_h(double x) {
            const double w = 1.0 / (x * x);
            return 1.7724538509055160273 * w * erfcx_tail_correction(w);
          }

          /** \ingroup opencl_kernels
           * h(y) - h(x) for x, y >= 4, with y - x = delta given exactly,
           * from the divided differences of the tail polynomials.
           *
           * @param x lower argument
           * @param y upper argument
           * @param delta y - x
           * @return h(y) - h(x)
           */
          double exp_mod_normal_erfcx_h_diff(double x, double y, double delta) {
            const double tp[6]
                = {0.000658749161529837803157, 0.0160837851487422766278,
                   0.125781726111229246204,    0.360344899949804439429,
                   0.305326634961232344035,    0.0163153871373020978498};
            const double tq[6]
                = {-0.00233520497626869185443, -0.0605183413124413191178,
                   -0.527905102951428412248,   -1.87295284992346047209,
                   -2.56852019228982242072,    -1.0};
            const double w0 = 1.0 / (x * x);
            const double w1 = 1.0 / (y * y);
            const double dw = -delta * (y + x) * w0 * w1;
            double a0 = 0.0;
            double q0 = 0.0;
            double q1 = 0.0;
            double dd_a = 0.0;
            double dd_q = 0.0;
            double pw0 = 1.0;
            double pw1 = 1.0;
            double sum_k = 0.0;
            for (int k = 0; k < 6; ++k) {
              q0 += tq[k] * pw0;
              q1 += tq[k] * pw1;
              if (k > 0) {
                dd_q += tq[k] * sum_k;
              }
              sum_k = sum_k * w0 + pw1;
              a0 += tp[k] * pw0 * w0;
              dd_a += tp[k] * sum_k;
              pw0 *= w0;
              pw1 *= w1;
            }
            return 1.7724538509055160273 * dw * (dd_a * q0 - a0 * dd_q)
                   / (q1 * q0);
          }

          /** \ingroup opencl_kernels
           * erfcx(lo) - erfcx(hi) for hi - lo = d > 0 with d given exactly,
           * and its log.
           *
           * @param lo lower argument
           * @param hi upper argument
           * @param d hi - lo
           * @param e_lo erfcx(lo)
           * @param e_hi erfcx(hi)
           * @param upper_exact nonzero if hi is the exact end
           * @param[out] log_delta log of the result
           * @return erfcx(lo) - erfcx(hi)
           */
          double exp_mod_normal_erfcx_diff(double lo, double hi, double d,
                                           double e_lo, double e_hi,
                                           int upper_exact, double* log_delta) {
            if (e_hi < 0.99 * e_lo) {
              const double delta = e_lo - e_hi;
              *log_delta = log(delta);
              return delta;
            }
            if (lo >= 4.0) {
              const double log_ratio_x
                  = upper_exact ? log1p(-d / hi) : -log1p(d / lo);
              const double log_r
                  = log_ratio_x
                    + log1p(exp_mod_normal_erfcx_h_diff(lo, hi, d)
                            / (1.0 + exp_mod_normal_erfcx_h(lo)));
              const double one_minus_r = -expm1(log_r);
              *log_delta = log(e_lo) + log(one_minus_r);
              return e_lo * one_minus_r;
            }
            const double m = upper_exact ? hi - 0.5 * d : lo + 0.5 * d;
            double e_km1 = erfcx(m);
            double e_k = erfcx_derivative(m, e_km1);
            double odd[4];
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              const double e_kp1 = 2.0 * m * e_k + 2.0 * k * e_km1;
              e_km1 = e_k;
              e_k = e_kp1;
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const double h = 0.5 * d;
            const double h2 = h * h;
            const double delta
                = -2.0 * h
                  * (odd[0]
                     + h2
                           * (odd[1] / 6.0
                              + h2 * (odd[2] / 120.0 + h2 * odd[3] / 5040.0)));
            *log_delta = log(delta);
            return delta;
          }

          /** \ingroup opencl_kernels
           * log(1 - x) for 0 <= x < 1.
           *
           * @param x argument
           * @return log(1 - x)
           */
          double exp_mod_normal_log1m(double x) {
            return x < 1e-17 ? -x : log1p(-x);
          }

          /** \ingroup opencl_kernels
           * The output `which` from F = g Delta, Delta = E(u) - E(t),
           * t >= 0, F <= 1/2.
           *
           * @return value (which = 0) or partial in z, v, sigma d/dsigma
           */
          double exp_mod_normal_lcdf_from_delta(
              double u, double t, double v, double log_g, double e_t,
              double delta, double log_delta, int upper, int which) {
            double value;
            double scale;
            if (upper) {
              const double g = exp(log_g);
              const double f = g * delta;
              value = exp_mod_normal_log1m(f);
              scale = -g / (1.0 - f);
            } else {
              value = log_g + log_delta;
              scale = delta > 1e-300 ? 1.0 / delta : exp(-log_delta);
            }
            if (which == 0) {
              return value;
            }
            const double de_t = erfcx_derivative(t, e_t);
            if (which == 1) {
              return scale * v * e_t;
            }
            if (which == 2) {
              return -scale * de_t * M_SQRT1_2;
            }
            return scale * v * M_SQRT2 * (u * e_t - 0.5 * de_t);
          }

          /** \ingroup opencl_kernels
           * The log ccdf where P = Phi(-z) + b is near underflow.
           *
           * @return value (which = 0) or partial in z, v, sigma d/dsigma
           */
          double exp_mod_normal_lccdf_tiny(
              double z, double v, double u, double t, double log_g,
              double log_b, double e_t, double e_mu, double e_mt, int which) {
            const double log_b_t
                = t >= 0.0 ? log_g + log(e_t)
                           : log_b + log1p(-0.5 * exp(-t * t) * e_mt);
            const double log_q = log_g + log(e_mu);
            const double log_p
                = fmax(log_q, log_b_t) + log1p(exp(-fabs(log_q - log_b_t)));
            if (which == 0) {
              return log_p;
            }
            const double b_p = exp(log_b_t - log_p);
            if (which == 1) {
              return -v * b_p;
            }
            const double g_p = exp(log_g - log_p);
            if (t >= 0.0) {
              const double de_t = erfcx_derivative(t, e_t);
              if (which == 2) {
                return g_p * de_t * M_SQRT1_2;
              }
              return -v * g_p * M_SQRT2 * (u * e_t - 0.5 * de_t);
            }
            const double phi_p = g_p * M_SQRT2 * M_2_SQRTPI * 0.5;
            if (which == 2) {
              return -(phi_p - (v - z) * b_p);
            }
            return -v * (phi_p - v * b_p);
          }

          /** \ingroup opencl_kernels
           * The log cdf (upper = 0) or log ccdf (upper = 1) of the
           * exponentially modified normal at z = (y - mu) / sigma and
           * v = lambda sigma, or one of its partials with respect to z, v
           * and sigma d/dsigma at fixed y, mu and lambda.
           *
           * @param z standardized argument
           * @param v lambda sigma
           * @param upper nonzero for the log ccdf
           * @param which 0: value, 1: d/dz, 2: d/dv, 3: sigma d/dsigma
           * @return the selected output
           */
          double exp_mod_normal_lcdf_part(double z, double v, int upper,
                                          int which) {
            const double u = -z * M_SQRT1_2;
            const double d = v * M_SQRT1_2;
            const double t = u + d;
            const double log_g = -M_LN2 - (0.5 * z) * z;
            double log_delta;
            if (z <= 0.0) {
              const double e_t = erfcx(t);
              const double delta = exp_mod_normal_erfcx_diff(
                  u, t, d, erfcx(u), e_t, 0, &log_delta);
              return exp_mod_normal_lcdf_from_delta(u, t, v, log_g, e_t, delta,
                                                    log_delta, upper, which);
            }
            const double log_b = d * (2.0 * u + d);
            if (t < -6.5 && log_b < -0.2877) {
              if (upper) {
                return which == 0   ? log_b
                       : which == 1 ? -v
                       : which == 2 ? v - z
                                    : v * v;
              }
              const double b = exp(log_b);
              if (which == 0) {
                return exp_mod_normal_log1m(b);
              }
              const double r = b / (1.0 - b);
              return which == 1 ? v * r : which == 2 ? (z - v) * r : -v * v * r;
            }
            const double g = exp(log_g);
            const double e_mu = erfcx(-u);
            const double q = g * e_mu;
            double e_t = 0.0;
            double e_mt = 0.0;
            double b;
            if (t >= 0.0) {
              e_t = erfcx(t);
              b = g * e_t;
            } else {
              e_mt = erfcx(-t);
              b = exp(log_b) - g * e_mt;
            }
            const double p = q + b;
            const int small_f = p > 0.75;
            if (small_f && t >= 0.0) {
              const double delta = exp_mod_normal_erfcx_diff(
                  u, t, d, erfcx(u), e_t, 0, &log_delta);
              return exp_mod_normal_lcdf_from_delta(u, t, v, log_g, e_t, delta,
                                                    log_delta, upper, which);
            }
            // the partial of F selected by which
            double f_x = 0.0;
            if (which == 1) {
              f_x = v * b;
            } else if (which > 1) {
              if (t >= 0.0) {
                const double de_t = erfcx_derivative(t, e_t);
                f_x = which == 2 ? -g * de_t * M_SQRT1_2
                                 : v * g * M_SQRT2 * (u * e_t - 0.5 * de_t);
              } else {
                const double phi = g * M_SQRT2 * M_2_SQRTPI * 0.5;
                f_x = which == 2 ? phi - (v - z) * b : v * (phi - v * b);
              }
            }
            if (small_f) {
              const double delta = exp_mod_normal_erfcx_diff(
                  -t, -u, d, e_mt, e_mu, 1, &log_delta);
              const double f = -expm1(log_b) + g * delta;
              if (!upper) {
                return which == 0 ? log(f) : f_x / f;
              }
              return which == 0 ? log1p(-f) : -f_x / (1.0 - f);
            }
            if (!upper) {
              const double f = 1.0 - p;
              if (which == 0) {
                return p < 0.25 ? exp_mod_normal_log1m(p) : log(f);
              }
              return f_x / f;
            }
            if (p > 1e-300) {
              return which == 0 ? log(p) : -f_x / p;
            }
            return exp_mod_normal_lccdf_tiny(z, v, u, t, log_g, log_b, e_t,
                                             e_mu, e_mt, which);
          }

          /** \ingroup opencl_kernels
           * Value and partials of the log cdf, one output per function.
           */
          double exp_mod_normal_lcdf_value(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 0, 0);
          } double exp_mod_normal_lcdf_d_z(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 0, 1);
          } double exp_mod_normal_lcdf_d_v(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 0, 2);
          } double exp_mod_normal_lcdf_d_s(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 0, 3);
          }

          /** \ingroup opencl_kernels
           * Value and partials of the log ccdf, one output per function.
           */
          double exp_mod_normal_lccdf_value(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 1, 0);
          } double exp_mod_normal_lccdf_d_z(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 1, 1);
          } double exp_mod_normal_lccdf_d_v(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 1, 2);
          } double exp_mod_normal_lccdf_d_s(double z, double v) {
            return exp_mod_normal_lcdf_part(z, v, 1, 3);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
