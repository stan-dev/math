#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_LOGIT_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_LOGIT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* inv_logit_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_LOGIT\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_INV_LOGIT\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           * Returns the inverse logit function applied to the kernel generator
           * expression.
           *
           * The inverse logit function is defined by
           *
           * \f$\mbox{logit}^{-1}(x) = \frac{1}{1 + \exp(-x)}\f$.
           *
           * This function can be used to implement the inverse link function
           * for logistic regression.
           *
           * The inverse to this function is <code>logit</code>.
           *
           * \f[
           * \mbox{inv\_logit}(y) =
           * \begin{cases}
           * \mbox{logit}^{-1}(y) & \mbox{if } -\infty\leq y \leq \infty \\[6pt]
           * \textrm{NaN} & \mbox{if } y = \textrm{NaN}
           * \end{cases}
           * \f]
           * \f[
           * \frac{\partial\, \mbox{inv\_logit}(y)}{\partial y} =
           * \begin{cases}
           * \frac{\partial\, \mbox{logit}^{-1}(y)}{\partial y} & \mbox{if }
           * -\infty\leq y\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } y =
           * \textrm{NaN} \end{cases} \f]
           *
           * \f[
           * \mbox{logit}^{-1}(y) = \frac{1}{1 + \exp(-y)}
           * \f]
           *
           * \f[
           * \frac{\partial \, \mbox{logit}^{-1}(y)}{\partial y} =
           * \frac{\exp(y)}{(\exp(y)+1)^2} \f]
           *
           * @param[in] x Argument.
           * @return natural logarithm of one minus the exponential of the
           * argument.
           */
          double inv_logit(double x) {
            if (x < 0) {
              if (x < log(2.2204460492503131E-16)) {
                return exp(x);
              }
              return exp(x) / (1 + exp(x));
            }
            return 1.0 / (1 + exp(-x));
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
