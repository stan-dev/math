#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static constexpr const char* beta_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BETA\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Return the beta function applied to the specified
           * arguments.
           *
           * For nonnegative arguments the value is exp(stan_lbeta(a, b)),
           * which keeps its accuracy for large arguments.  Negative and
           * NaN arguments use the lgamma difference.  The kernel must
           * include the lgamma_stirling, lgamma_stirling_diff and lbeta
           * device functions before this one.
           *
           * @param a First value
           * @param b Second value
           * @return beta(a, b)
           */
          double beta(double a, double b) {
            if (a >= 0 && b >= 0) {
              return exp(stan_lbeta(a, b));
            }
            return exp(lgamma(a) + lgamma(b) - lgamma(a + b));
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
