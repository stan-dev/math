#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGIT_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGIT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* logit_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGIT\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOGIT\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_kernels
           *
           * Return the log odds applied to the kernel generator
           * expression.
           *
           * The logit function is defined as for \f$x \in [0, 1]\f$ by
           * returning the log odds of \f$x\f$ treated as a probability,
           *
           * \f$\mbox{logit}(x) = \log \left( \frac{x}{1 - x} \right)\f$.
           *
           * The inverse to this function is <code>inv_logit</code>.
           *
           *
           \f[
           \mbox{logit}(x) =
          \begin{cases}
          \textrm{NaN}& \mbox{if } x < 0 \textrm{ or } x > 1\\
          \ln\frac{x}{1-x} & \mbox{if } 0\leq x \leq 1 \\[6pt]
          \textrm{NaN} & \mbox{if } x = \textrm{NaN}
          \end{cases}
          \f]

          \f[
          \frac{\partial\, \mbox{logit}(x)}{\partial x} =
          \begin{cases}
          \textrm{NaN}& \mbox{if } x < 0 \textrm{ or } x > 1\\
          \frac{1}{x-x^2}& \mbox{if } 0\leq x\leq 1 \\[6pt]
          \textrm{NaN} & \mbox{if } x = \textrm{NaN}
          \end{cases}
          \f]
          *
          * @param x argument
          * @return log odds of argument
          */
          double logit(double x) { return log(x) - log1m(x); }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
