#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LBETA_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LBETA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* lbeta_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LBETA\n"
      "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LBETA\n" STRINGIFY(
          // \endcond
          /**
           * Return the log of the beta function applied to the specified
           * arguments.
           *
           * The beta function is defined for \f$a > 0\f$ and \f$b > 0\f$ by
           *
           * \f$\mbox{B}(a, b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}\f$.
           *
           * This function returns its log,
           *
           * \f$\log \mbox{B}(a, b) = \log \Gamma(a) + \log \Gamma(b) - \log
           \Gamma(a+b)\f$.
           *
           * See stan::math::lgamma() for the double-based and stan::math for
           the
           * variable-based log Gamma function.
           * This function is numerically more stable than naive evaluation via
           lgamma.
           *
             \f[
             \mbox{lbeta}(\alpha, \beta) =
             \begin{cases}
               \ln\int_0^1 u^{\alpha - 1} (1 - u)^{\beta - 1} \, du & \mbox{if }
           \alpha, \beta>0 \\[6pt] \textrm{NaN} & \mbox{if } \alpha =
           \textrm{NaN or } \beta = \textrm{NaN} \end{cases} \f]

             \f[
             \frac{\partial\, \mbox{lbeta}(\alpha, \beta)}{\partial \alpha} =
             \begin{cases}
               \Psi(\alpha)-\Psi(\alpha+\beta) & \mbox{if } \alpha, \beta>0
           \\[6pt] \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta =
           \textrm{NaN} \end{cases} \f]

             \f[
             \frac{\partial\, \mbox{lbeta}(\alpha, \beta)}{\partial \beta} =
             \begin{cases}
               \Psi(\beta)-\Psi(\alpha+\beta) & \mbox{if } \alpha, \beta>0
           \\[6pt] \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta =
           \textrm{NaN} \end{cases} \f]
           *
           * @param a First value
           * @param b Second value
           * @return Log of the beta function applied to the two values.
           */
          double lbeta(double a, double b) {
            if (isnan(a) || isnan(b)) {
              return a;
            }

            double x;  // x is the smaller of the two
            double y;
            if (a < b) {
              x = a;
              y = b;
            } else {
              x = b;
              y = a;
            }

            // Special cases
            if (x == 0) {
              return INFINITY;
            }
            if (isinf(y)) {
              return -INFINITY;
            }

            // For large x or y, separate the lgamma values into Stirling
            // approximations and appropriate corrections. The Stirling
            // approximations allow for analytic simplification and the
            // corrections are added later.
            //
            // The overall approach is inspired by the code in R, where the
            // algorithm is credited to W. Fullerton of Los Alamos Scientific
            // Laboratory
            if (y < LGAMMA_STIRLING_DIFF_USEFUL) {
              // both small
              return lgamma(x) + lgamma(y) - lgamma(x + y);
            }
            double x_over_xy = x / (x + y);
            if (x < LGAMMA_STIRLING_DIFF_USEFUL) {
              // y large, x small
              double stirling_diff
                  = lgamma_stirling_diff(y) - lgamma_stirling_diff(x + y);
              double stirling
                  = (y - 0.5) * log1p(-x_over_xy) + x * (1 - log(x + y));
              return stirling + lgamma(x) + stirling_diff;
            }

            // both large
            double stirling_diff = lgamma_stirling_diff(x)
                                   + lgamma_stirling_diff(y)
                                   - lgamma_stirling_diff(x + y);
            double stirling = (x - 0.5) * log(x_over_xy) + y * log1p(-x_over_xy)
                              + 0.5 * log(2.0 * M_PI) - 0.5 * log(y);
            return stirling + stirling_diff;
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
