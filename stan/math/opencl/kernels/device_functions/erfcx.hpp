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
           * Return the scaled complementary error function
           * exp(x * x) * erfc(x) of the kernel generator expression.
           *
           * Mirrors stan/math/prim/fun/erfcx.hpp: the Cody (1969) rational
           * approximation above 4, and below it exp(x * x) * erfc(x) with
           * x * x Dekker-split so that exp() receives an exact argument.
           * Keep the two in step; the split is what holds the lower branch
           * to a few ulp instead of x * x * eps.
           *
           * @param x argument
           * @return scaled complementary error function of the argument
           */
          double erfcx(double x) {
            if (x >= 4.0) {
              double u = 1.0 / (x * x);
              double p
                  = 0.000658749161529837803157
                    + u
                          * (0.0160837851487422766278
                             + u
                                   * (0.125781726111229246204
                                      + u
                                            * (0.360344899949804439429
                                               + u
                                                     * (0.305326634961232344035
                                                        + u * 0.0163153871373020978498))));
              double q
                  = -0.00233520497626869185443
                    + u
                          * (-0.0605183413124413191178
                             + u
                                   * (-0.527905102951428412248
                                      + u
                                            * (-1.87295284992346047209
                                               + u
                                                     * (-2.56852019228982242072
                                                        + u * -1.0))));
              return (M_2_SQRTPI * 0.5 + (p / q) * u) / x;
            }
            if (x < -27.0) {
              return INFINITY;
            }
            double t = 134217729.0 * x;
            double x_hi = t - (t - x);
            double x_lo = x - x_hi;
            return exp(x_hi * x_hi) * exp(x_lo * (x + x_hi)) * erfc(x);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
