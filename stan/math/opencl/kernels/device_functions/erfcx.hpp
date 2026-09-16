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
           * the rounding of x * x corrected so that exp() effectively
           * receives an exact argument. Without that correction the error
           * of x * x is amplified by exp into roughly x * x * eps, which is
           * 512 ulp at x = -26.
           *
           * The correction is written with fma() rather than the Dekker
           * split used on the host. The split relies on t - (t - x) not
           * being simplified to x, which is true in floating point but not
           * over the reals, and the OpenCL compiler simplifies it: measured
           * on a Tesla V100 (OpenCL 3.0 CUDA), the split form returned
           * x_lo == 0 at every one of 4096 test points and scored 512.66
           * ulp, bit-identical to the uncorrected formula. fma() is a
           * single instruction, so there is nothing to reassociate; it
           * measures 3.40 ulp on the same device.
           *
           * @param x argument
           * @return scaled complementary error function of the argument
           */
          double erfcx(double x) {
            if (x >= 4.0) {
              double u = 1.0 / (x * x);
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
              return (M_2_SQRTPI * 0.5 + (p / q) * u) / x;
            }
            if (x < -27.0) {
              return INFINITY;
            }
            double h = x * x;
            return exp(h) * (1.0 + fma(x, x, -h)) * erfc(x);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
