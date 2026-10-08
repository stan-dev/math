#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA_DIFF_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA_DIFF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* digamma_diff_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA_DIFF\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA_DIFF\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_device_functions
           * Return the difference of the digamma function at two arguments
           * that differ by a nonnegative offset, psi(x + d) - psi(x).
           *
           * The plain difference of two digamma calls loses all accuracy
           * when x is large. This function has a relative error of a few
           * ulp for all x > 0 and d >= 0. It uses the same method as
           * stan::math::digamma_diff(): for an integer d from 0 to 8 the
           * sum of 1 / (x + j), j < d; for x < 10 and d >= 10 the plain
           * difference, which does not cancel much there; otherwise shift x
           * up to y >= 10 with the recurrence psi(y + 1) = psi(y) + 1 / y,
           * then use the asymptotic expansion of psi(y + d) - psi(y), with
           * the differences of the powers of 1 / y^2 and 1 / (y + d)^2
           * formed without cancellation. Needs the digamma device function.
           *
           * @param x first argument, positive
           * @param d offset, nonnegative
           * @return psi(x + d) - psi(x), or NaN if x is not positive or d
           * is negative
           */
          double digamma_diff(double x, double d) {
            if (isnan(x) || isnan(d) || !(x > 0) || d < 0) {
              return NAN;
            }
            if (isinf(d)) {
              return INFINITY;
            }
            if (isinf(x)) {
              return 0.0;
            }
            // a count d from 0 to 8: sum of the positive terms 1 / (x + j),
            // from the smallest (exactly 1 / x for d = 1)
            if (d <= 8.0 && d == floor(d)) {
              double sum = 0.0;
              for (int j = (int)d - 1; j >= 0; --j) {
                sum += 1.0 / (x + j);
              }
              return sum;
            }
            // x < 10 and d >= 10: the plain difference does not cancel much
            if (x < 10.0 && d >= 10.0) {
              return digamma(x + d) - digamma(x);
            }
            // B_{2i} / (2i), i = 1..8
            const double coeffs[8] = {
                1.0 / 12.0,  -1.0 / 120.0,     1.0 / 252.0, -1.0 / 240.0,
                1.0 / 132.0, -691.0 / 32760.0, 1.0 / 12.0,  -3617.0 / 8160.0};
            // each shift term is d / (y (y + d)) = 1 / y - 1 / (y + d); for
            // d < y it is formed as (d / (y + d)) / y, which neither cancels
            // nor overflows; for d >= y the parts 1 / y are summed separately
            // and added last, so that for small x the dominant 1 / x keeps its
            // correct rounding (for d = 1 the result is exactly 1 / x)
            double inv_sum = 0.0;
            double shift_sum = 0.0;
            double y = x;
            while (y < 10.0) {
              if (d >= y) {
                inv_sum += 1.0 / y;
                shift_sum -= 1.0 / (y + d);
              } else {
                shift_sum += (d / (y + d)) / y;
              }
              y += 1.0;
            }
            const double y_plus_d = y + d;
            const double d_frac = d / y_plus_d;
            // the square of the inverse underflows to 0 where the inverse of
            // the square would overflow
            const double inv_y = 1.0 / y;
            const double inv_y_plus_d = 1.0 / y_plus_d;
            const double u = inv_y * inv_y;
            const double v = inv_y_plus_d * inv_y_plus_d;
            // u - v = u (d / (y + d)) (1 + y / (y + d)), without cancellation
            const double u_minus_v = u * d_frac * (1.0 + y / y_plus_d);
            // u^i - v^i = (u - v) h_i with h_1 = 1, h_(i+1) = u h_i + v^i
            double h = 1.0;
            double v_pow = v;
            double series = coeffs[0];
            for (int i = 1; i < 8; ++i) {
              h = u * h + v_pow;
              v_pow *= v;
              series += coeffs[i] * h;
            }
            return inv_sum
                   + (shift_sum + log1p(d / y) + 0.5 * d_frac / y
                      + u_minus_v * series);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
