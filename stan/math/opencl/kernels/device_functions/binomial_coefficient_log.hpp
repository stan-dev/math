#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BINOMIAL_COEFFICIENT_LOG_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BINOMIAL_COEFFICIENT_LOG_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* binomial_coefficient_log_device_function
    = "\n"
      "#ifndef "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BINOMIAL_COEFFICIENT_LOG\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_BINOMIAL_COEFFICIENT_"
      "LOG\n" STRINGIFY(
          // \endcond
          /**
           * Return the log of the binomial coefficient for the specified
           * arguments.
           *
           * The binomial coefficient, \f${n \choose k}\f$, read "n choose k",
           is
           * defined for \f$0 \leq k \leq n\f$ by
           *
           * \f${n \choose k} = \frac{n!}{k! (n-k)!}\f$.
           *
           * This function uses Gamma functions to define the log
           * and generalize the arguments to continuous n and k.
           *
           * \f$ \log {n \choose k}
           * = \log \ \Gamma(n+1) - \log \Gamma(k+1) - \log \Gamma(n-k+1)\f$.
           *
           *
             \f[
             \mbox{binomial\_coefficient\_log}(x, y) =
             \begin{cases}
               \textrm{error} & \mbox{if } y > x + 1 \textrm{ or } y < -1
           \textrm{ or } x
           < -1\\
               \ln\Gamma(x+1) & \mbox{if } -1 < y < x + 1 \\
               \quad -\ln\Gamma(y+1)& \\
               \quad -\ln\Gamma(x-y+1)& \\[6pt]
               \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
             \end{cases}
             \f]

             \f[
             \frac{\partial\, \mbox{binomial\_coefficient\_log}(x, y)}{\partial
           x} = \begin{cases} \textrm{error} & \mbox{if } y > x + 1 \textrm{ or
           } y < -1 \textrm{ or } x
           < -1\\
               \Psi(x+1) & \mbox{if } 0\leq y \leq x \\
               \quad -\Psi(x-y+1)& \\[6pt]
               \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
             \end{cases}
             \f]

             \f[
             \frac{\partial\, \mbox{binomial\_coefficient\_log}(x, y)}{\partial
           y} = \begin{cases} \textrm{error} & \mbox{if } y > x + 1 \textrm{ or
           } y < -1 \textrm{ or } x
           < -1\\
               -\Psi(y+1) & \mbox{if } 0\leq y \leq x \\
               \quad +\Psi(x-y+1)& \\[6pt]
               \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
             \end{cases}
             \f]
           *
           *  This function is numerically more stable than naive evaluation via
           lgamma.
           *
           * @param n total number of objects.
           * @param k number of objects chosen.
           * @return log (n choose k).
           */
          double binomial_coefficient_log(double n, double k) {
            if (isnan(n) || isnan(k)) {
              return NAN;
            }

            // Choosing the more stable of the symmetric branches
            if (n > -1 && k > n / 2.0 + 1e-8) {
              k = n - k;
            }

            double n_plus_1 = n + 1;
            double n_plus_1_mk = n_plus_1 - k;

            if (k == 0) {
              return 0;
            } else if (n_plus_1 < LGAMMA_STIRLING_DIFF_USEFUL) {
              return lgamma(n_plus_1) - lgamma(k + 1) - lgamma(n_plus_1_mk);
            } else {
              return -lbeta(n_plus_1_mk, k + 1) - log1p(n);
            }
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif  // BINOMIAL_COEFFICIENT_LOG_HPP
