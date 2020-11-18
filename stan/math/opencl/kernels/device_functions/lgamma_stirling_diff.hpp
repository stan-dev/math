#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING_DIFF_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING_DIFF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* lgamma_stirling_diff_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING_DIFF\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LGAMMA_STIRLING_DIFF\n"
      "#define LGAMMA_STIRLING_DIFF_USEFUL 10\n" STRINGIFY(
          // \endcond
          /**
           * Return the difference between log of the gamma function and its
           Stirling
           * approximation.
           * This is useful to stably compute log of ratios of gamma functions
           with large
           * arguments where the Stirling approximation allows for analytic
           solution
           * and the (small) differences can be added afterwards.
           * This is for example used in the implementation of lbeta.
           *
           * The function will return correct value for all arguments, but using
           it can
           * lead to a loss of precision when x < lgamma_stirling_diff_useful.
           *
             \f[
             \mbox{lgamma_stirling_diff}(x) =
              \log(\Gamma(x)) - \frac{1}{2} \log(2\pi) +
              (x-\frac{1}{2})*\log(x) - x
             \f]

           *
           * @param x value
           * @return Difference between lgamma(x) and its Stirling
           approximation.
           */
          double lgamma_stirling_diff(double x) {
            if (isnan(x)) {
              return x;
            }
            if (x == 0) {
              return INFINITY;
            }
            if (x < LGAMMA_STIRLING_DIFF_USEFUL) {
              return lgamma(x) - lgamma_stirling(x);
            }
            // Using the Stirling series as expressed in formula 5.11.1. at
            // https://dlmf.nist.gov/5.11
            const double stirling_series[] = {
                0.0833333333333333333333333,   -0.00277777777777777777777778,
                0.000793650793650793650793651, -0.000595238095238095238095238,
                0.000841750841750841750841751, -0.00191752691752691752691753,
                0.00641025641025641025641026,  -0.0295506535947712418300654};
            double multiplier = 1 / x;
            double inv_x_squared = multiplier * multiplier;
            double result = stirling_series[0] * multiplier;
            for (int n = 1; n < 6; n++) {
              multiplier *= inv_x_squared;
              result += stirling_series[n] * multiplier;
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
