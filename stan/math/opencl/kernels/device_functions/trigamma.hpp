#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_TRIGAMMA_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_TRIGAMMA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* trigamma_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_TRIGAMMA\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_TRIGAMMA\n" STRINGIFY(
          // \endcond
          /**
           * Return the trigamma function applied to the argument.   The
           * trigamma function returns the second derivative of the
           * log Gamma function evaluated at the specified argument.
           *
           * @param x argument
           * @return second derivative of log Gamma function at argument
           *
           */
          double trigamma(double x) {
            // negative integers and zero return positive infinity
            // see http://mathworld.wolfram.com/PolygammaFunction.html
            if (x <= 0.0 && floor(x) == x) {
              return INFINITY;
            }

            // negative non-integers: use the reflection formula
            // see http://mathworld.wolfram.com/PolygammaFunction.html
            double value = 0;
            bool neg = false;
            if (x <= 0 && floor(x) != x) {
              neg = true;
              double v_sqrt = M_PI / sinpi(-x);
              value = -v_sqrt * v_sqrt;
              x = -x + 1.0;
            }

            double small = 0.0001;
            double large = 5.0;
            // small value approximation if x <= small.
            if (x <= small) {
              value = 1.0 / (x * x);
            } else {
              // use recurrence relation until x >= large
              // see http://mathworld.wolfram.com/PolygammaFunction.html
              double z = x;
              while (z < large) {
                value += 1.0 / (z * z);
                z += 1.0;
              }

              // bernoulli numbers
              double b2 = 1.0 / 6.0;
              double b4 = -1.0 / 30.0;
              double b6 = 1.0 / 42.0;
              double b8 = -1.0 / 30.0;

              // asymptotic expansion as a Laurent series if x >= large
              // see http://en.wikipedia.org/wiki/Trigamma_function
              double y = 1.0 / (z * z);
              value += 0.5 * y
                       + (1.0 + y * (b2 + y * (b4 + y * (b6 + y * b8)))) / z;
            }
            if (neg) {
              return -value;
            } else {
              return value;
            }
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
