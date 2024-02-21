#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_INV_LOGIT_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_INV_LOGIT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* log1m_inv_logit_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_INV_LOGIT\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG1M_INV_LOGIT\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           *
           * Return the the natural logarithm of 1 minus the inverse logit
           * applied to the kernel generator expression.
           *
             \f[
             \mbox{log1m\_inv\_logit}(x) =
            \begin{cases}
              -\ln(\exp(x)+1) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
              \textrm{NaN} & \mbox{if } x = \textrm{NaN}
            \end{cases}
            \f]

            \f[
            \frac{\partial\, \mbox{log1m\_inv\_logit}(x)}{\partial x} =
            \begin{cases}
              -\frac{\exp(x)}{\exp(x)+1} & \mbox{if } -\infty\leq x\leq \infty
          \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN} \end{cases} \f]
          *
          * @param x argument
          * @return log of one minus the inverse logit of the argument
          */
          inline double log1m_inv_logit(double x) {
            if (x > 0.0) {
              return -x - log1p_exp(-x);  // prevent underflow
            }
            return -log1p_exp(x);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
