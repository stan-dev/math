#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGISTIC_TAIL_DERIV_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGISTIC_TAIL_DERIV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* logistic_tail_deriv_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGISTIC_TAIL_DERIV\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGISTIC_TAIL_"
      "DERIV\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           *
           * Return inv_logit(-x) / sigma.
           *
           * @param x scaled difference
           * @param sigma scale
           * @return inv_logit(-x) / sigma
           */
          double logistic_tail_deriv(double x, double sigma) {
            if (x > 700.0) {
              return exp(log1m_inv_logit(x) - log(sigma));
            }
            return inv_logit(-x) / sigma;
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
