#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_ERFCX_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_ERFCX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* erfcx_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_ERFCX\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_ERFCX\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           *
           * Correction factor of the Cody (1969) third-interval rational:
           * erfcx(x) = (INV_SQRT_PI + u * correction(u)) / x, u = 1 / x^2.
           *
           * Split out so the derivative can reuse it. See
           * erfcx_tail_derivative.
           *
           * @param u inverse square of the argument, 0 <= u <= 1/16
           * @return P(u) / Q(u)
           */
          double erfcx_tail_correction(double u) {
            double p = 0.0163153871373020978498;
            p = 0.305326634961232344035 + u * p;
            p = 0.360344899949804439429 + u * p;
            p = 0.125781726111229246204 + u * p;
            p = 0.0160837851487422766278 + u * p;
            p = 0.000658749161529837803157 + u * p;
            double q = -1.0;
            q = -2.56852019228982242072 + u * q;
            q = -1.87295284992346047209 + u * q;
            q = -0.527905102951428412248 + u * q;
            q = -0.0605183413124413191178 + u * q;
            q = -0.00233520497626869185443 + u * q;
            return p / q;
          }

          /** \ingroup opencl_kernels
           *
           * Cody (1969) third-interval rational, for x >= 4.
           *
           * @param x argument
           * @return scaled complementary error function
           */
          double erfcx_cody_tail(double x) {
            double u = 1.0 / (x * x);
            return (M_2_SQRTPI * 0.5 + u * erfcx_tail_correction(u)) / x;
          }

          /** \ingroup opencl_kernels
           *
           * Derivative of erfcx, 2 * x * erfcx(x) - 2 / sqrt(pi).
           *
           * That difference cancels for large x: both terms approach
           * 2 / sqrt(pi) while the result decays like 1 / (sqrt(pi) * x^2).
           * The difference form reaches 2.55e+11 ulp at x = 1e6. For
           * x >= 4 the constant cancels analytically against the leading
           * term of the tail rational, leaving 2 * u * C(u).
           *
           * Takes the value as an argument so the reverse pass does not
           * evaluate erfcx twice.
           *
           * @param x argument
           * @param value erfcx(x)
           * @return derivative of erfcx at x
           */
          double erfcx_derivative(double x, double value) {
            if (x >= 4.0) {
              double u = 1.0 / (x * x);
              return 2.0 * u * erfcx_tail_correction(u);
            }
            return 2.0 * x * value - M_2_SQRTPI;
          }

          /** \ingroup opencl_kernels
           *
           * Cody (1969) second-interval rational, for 0.46875 <= x <= 4.
           * Yields erfcx directly: the exponential is cancelled
           * analytically, so there is no exp and no erfc call.
           *
           * @param y argument
           * @return scaled complementary error function
           */
          double erfcx_cody_middle(double y) {
            double p = 2.15311535474403846e-8 * y;
            p = (p + 5.64188496988670089e-1) * y;
            p = (p + 8.88314979438837594) * y;
            p = (p + 66.1191906371416295) * y;
            p = (p + 298.635138197400131) * y;
            p = (p + 881.952221241769090) * y;
            p = (p + 1712.04761263407058) * y;
            p = (p + 2051.07837782607147) * y;
            double q = y;
            q = (q + 15.7449261107098347) * y;
            q = (q + 117.693950891312499) * y;
            q = (q + 537.181101862009858) * y;
            q = (q + 1621.38957456669019) * y;
            q = (q + 3290.79923573345963) * y;
            q = (q + 4362.61909014324716) * y;
            q = (q + 3439.36767414372164) * y;
            return (p + 1230.33935479799725) / (q + 1230.33935480374942);
          }

          /** \ingroup opencl_kernels
           *
           * Degree-18 Chebyshev-economized expansion of erfcx about zero,
           * for |x| < 0.46875. Covers both signs with no branch and no
           * library call.
           *
           * @param x argument
           * @return scaled complementary error function
           */
          double erfcx_small(double x) {
            double p = 3.05977060678449757e-06;
            p = -9.35890030086883823e-06 + x * p;
            p = 2.46655529768908249e-05 + x * p;
            p = -7.08163358203131886e-05 + x * p;
            p = 1.98445679338826757e-04 + x * p;
            p = -5.34506929034156810e-04 + x * p;
            p = 1.38888415444527033e-03 + x * p;
            p = -3.47359067853470795e-03 + x * p;
            p = 8.33333374332981443e-03 + x * p;
            p = -1.91048337772546720e-02 + x * p;
            p = 4.16666666458337179e-02 + x * p;
            p = -8.59717459974174147e-02 + x * p;
            p = 1.66666666667239644e-01 + x * p;
            p = -3.00901111227312890e-01 + x * p;
            p = 4.99999999999992839e-01 + x * p;
            p = -7.52252778063651983e-01 + x * p;
            p = 1.0 + x * p;
            p = -1.12837916709551256 + x * p;
            return 1.0 + x * p;
          }

          /** \ingroup opencl_kernels
           *
           * Return the scaled complementary error function
           * exp(x * x) * erfc(x) of the kernel generator expression.
           *
           * Mirrors stan/math/prim/fun/erfcx.hpp branch for branch and
           * formula for formula; the two must be kept in step. The whole
           * positive axis is covered without a library call, which is
           * what makes this fast. On the negative side exp is
           * unavoidable, since erfcx grows like 2*exp(x*x); its argument
           * is corrected with fma, because exp amplifies the rounding of
           * x * x into roughly x * x * eps, which is 512 ulp at x = -26.
           *
           * The correction must be written with fma and not with a Dekker
           * split. The split relies on t - (t - x) not being simplified
           * to x, which holds in floating point but not over the reals,
           * and the OpenCL compiler does simplify it: the split form
           * silently degrades to the uncorrected product on device while
           * still being correct on the host. Do not reintroduce it.
           *
           * @param x argument
           * @return scaled complementary error function of the argument
           */
          double erfcx(double x) {
            if (x >= 4.0) {
              return erfcx_cody_tail(x);
            }
            if (x >= 0.46875) {
              return erfcx_cody_middle(x);
            }
            if (x > -0.46875) {
              return erfcx_small(x);
            }
            if (x < -27.0) {
              return INFINITY;
            }
            double h = x * x;
            double two_exp_x2 = 2.0 * exp(h) * (1.0 + fma(x, x, -h));
            if (x < -6.1) {
              return two_exp_x2;
            }
            double y = -x;
            return two_exp_x2
                   - (y >= 4.0 ? erfcx_cody_tail(y) : erfcx_cody_middle(y));
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
