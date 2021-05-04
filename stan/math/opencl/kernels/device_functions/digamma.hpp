#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const char* digamma_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_DIGAMMA\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Calculates the digamma function - derivative of logarithm of gamma.
           * This implementation is based on one from boost 1.69.0:
           * https://www.boost.org/doc/libs/1_69_0/boost/math/special_functions/digamma.hpp.
           *
           * @param x point at which to calculate digamma
           * @return digamma(x)
           */
          double digamma(double x) {
            double result = 0;
            if (x <= -1) {
              x = 1 - x;
              double remainder = x - floor(x);
              if (remainder > 0.5) {
                remainder -= 1;
              }
              if (remainder == 0) {
                return NAN;
              }
              result = M_PI / tan(M_PI * remainder);
            }
            if (x == 0) {
              return NAN;
            }
            // in boost: x >= digamma_large_lim(t)
            if (x > 10) {
              // in boost: result += digamma_imp_large(x, t);
              const double P[8]
                  = {0.083333333333333333333333333333333333333333333333333,
                     -0.0083333333333333333333333333333333333333333333333333,
                     0.003968253968253968253968253968253968253968253968254,
                     -0.0041666666666666666666666666666666666666666666666667,
                     0.0075757575757575757575757575757575757575757575757576,
                     -0.021092796092796092796092796092796092796092796092796,
                     0.083333333333333333333333333333333333333333333333333,
                     -0.44325980392156862745098039215686274509803921568627};
              x -= 1;
              result += log(x);
              result += 1 / (2 * x);
              double z = 1 / (x * x);
              double tmp = P[7];
              for (int i = 6; i >= 0; i--) {
                tmp = tmp * z + P[i];
              }
              // tmp=boost::tools::evaluate_polynomial(P, z);
              result -= z * tmp;
            } else {
              while (x > 2) {
                x -= 1;
                result += 1 / x;
              }
              while (x < 1) {
                result -= 1 / x;
                x += 1;
              }
              // in boost: result += digamma_imp_1_2(x, t);
              const float Y = 0.99558162689208984F;

              const double root1 = (double)1569415565 / 1073741824uL;  // NOLINT
              const double root2
                  = (double)381566830 / 1073741824uL / 1073741824uL;  // NOLINT
              const double root3
                  = 0.9016312093258695918615325266959189453125e-19;

              const double P[6]
                  = {0.25479851061131551,   -0.32555031186804491,
                     -0.65031853770896507,  -0.28919126444774784,
                     -0.045251321448739056, -0.0020713321167745952};
              const double Q[7] = {1.0,
                                   2.0767117023730469,
                                   1.4606242909763515,
                                   0.43593529692665969,
                                   0.054151797245674225,
                                   0.0021284987017821144,
                                   -0.55789841321675513e-6};
              double g = x - root1 - root2 - root3;
              double tmp = P[5];
              for (int i = 4; i >= 0; i--) {
                tmp = tmp * (x - 1) + P[i];
              }
              double tmp2 = Q[6];
              for (int i = 5; i >= 0; i--) {
                tmp2 = tmp2 * (x - 1) + Q[i];
              }
              // in boost: T r = tools::evaluate_polynomial(P, T(x-1)) /
              // tools::evaluate_polynomial(Q, T(x-1));
              double r = tmp / tmp2;
              result += g * Y + g * r;
            }
            return result;
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
