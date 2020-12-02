#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_MULTIPLY_LOG_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_MULTIPLY_LOG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* multiply_log_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_MULTIPLY_LOG\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_MULTIPLY_LOG\n" STRINGIFY(
          // \endcond
          /**
           * Calculate the value of the first argument
           * times log of the second argument while behaving
           * properly with 0 inputs.
           *
           * \f$ a * \log b \f$.
           *
             \f[
             \mbox{multiply\_log}(x, y) =
             \begin{cases}
               0 & \mbox{if } x=y=0\\
               x\ln y & \mbox{if } x, y\neq0 \\[6pt]
               \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
             \end{cases}
             \f]

             \f[
             \frac{\partial\, \mbox{multiply\_log}(x, y)}{\partial x} =
             \begin{cases}
               \infty & \mbox{if } x=y=0\\
               \ln y & \mbox{if } x, y\neq 0 \\[6pt]
               \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
             \end{cases}
             \f]

             \f[
             \frac{\partial\, \mbox{multiply\_log}(x, y)}{\partial y} =
             \begin{cases}
               \infty & \mbox{if } x=y=0\\
               \frac{x}{y} & \mbox{if } x, y\neq 0 \\[6pt]
               \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
             \end{cases}
             \f]
           *
           * @param a the first variable
           * @param b the second variable
           * @return a * log(b)
           */
          double multiply_log(double a, double b) {
            if (b == 0.0 && a == 0.0) {
              return 0.0;
            }
            return a * log(b);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
