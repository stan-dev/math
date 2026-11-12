#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_BETA_RATIO_HPP
#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_BETA_RATIO_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static constexpr const char* log_beta_ratio_device_function
    = "\n"
      "#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_BETA_RATIO\n"
      "#define "
      "STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_LOG_BETA_RATIO\n" STRINGIFY(
          // \endcond
          /** \ingroup opencl_device_functions
           * Return log1p(t) - t for 0 <= t < 1 without cancellation. See
           * stan::math::internal::log1pmx().
           *
           * @param t argument in [0, 1)
           * @return log1p(t) - t
           */
          double stan_log1pmx(double t) {
            if (t > 0.25) {
              return log1p(t) - t;
            }
            const double r = t / (2.0 + t);
            const double r2 = r * r;
            double s = 1.0 / 21.0;
            s = 1.0 / 19.0 + r2 * s;
            s = 1.0 / 17.0 + r2 * s;
            s = 1.0 / 15.0 + r2 * s;
            s = 1.0 / 13.0 + r2 * s;
            s = 1.0 / 11.0 + r2 * s;
            s = 1.0 / 9.0 + r2 * s;
            s = 1.0 / 7.0 + r2 * s;
            s = 1.0 / 5.0 + r2 * s;
            s = 1.0 / 3.0 + r2 * s;
            return r * (2.0 * r2 * s - t);
          }

          /** \ingroup opencl_device_functions
           * Return the count that stan_log_beta_ratio_term() removes from
           * its term: k if 0 < k < x, else 0.
           *
           * @param x shape
           * @param k count, a nonnegative integer
           * @return the removed count
           */
          double stan_log_beta_ratio_removed(double x, double k) {
            if (k == 0 || !(k / x < 1.0)) {
              return 0.0;
            }
            return k;
          }

          /** \ingroup opencl_device_functions
           * Return (x - 1/2) log1p(k / x) minus the count that
           * stan_log_beta_ratio_removed() returns. For k < x the term is
           * about k and is formed as x (log1p(k / x) - k / x)
           * - log1p(k / x) / 2. See
           * stan::math::internal::stirling_log1p_term().
           *
           * @param x shape, at least LGAMMA_STIRLING_DIFF_USEFUL
           * @param k count, a nonnegative integer
           * @return the term without the removed count
           */
          double stan_log_beta_ratio_term(double x, double k) {
            if (k == 0) {
              return 0.0;
            }
            const double t = k / x;
            if (t < 1.0) {
              return x * stan_log1pmx(t) - 0.5 * log1p(t);
            }
            return (x - 0.5) * log1p(t);
          }

          /** \ingroup opencl_device_functions
           * Return an estimate of the size of the terms that lbeta(a, b)
           * adds up: min(a, b) (1 + log1p(max / min)). Only used to choose
           * between two forms.
           *
           * @param a first argument
           * @param b second argument
           * @return the size estimate
           */
          double stan_lbeta_terms_size(double a, double b) {
            const double small = fmin(a, b);
            const double large = fmax(a, b);
            return small * (1.0 + log1p(large / small));
          }

          /** \ingroup opencl_device_functions
           * Return lbeta(alpha + n, beta + m) - lbeta(alpha, beta) for
           * shapes alpha, beta > 0 and integer counts n, m >= 0, without
           * the cancellation of the two lbeta values for large shapes.
           * This is
           * stan::math::internal::log_beta_ratio() with its shape-only
           * part (internal::log_beta_ratio_denominator()) computed inside.
           *
           * For shapes of at least LGAMMA_STIRLING_DIFF_USEFUL, every lgamma
           * is written as its Stirling form plus lgamma_stirling_diff().
           * The terms linear in the arguments cancel exactly. The function
           * takes this form or the lbeta form, whichever has the smaller
           * terms.
           *
           * @param alpha first shape, positive
           * @param beta second shape, positive
           * @param n first count, a nonnegative integer
           * @param m second count, a nonnegative integer
           * @return lbeta(alpha + n, beta + m) - lbeta(alpha, beta)
           */
          double stan_log_beta_ratio(double alpha, double beta, double n,
                                     double m) {
            if (!(alpha >= LGAMMA_STIRLING_DIFF_USEFUL
                  && beta >= LGAMMA_STIRLING_DIFF_USEFUL)) {
              return stan_lbeta(alpha + n, beta + m) - stan_lbeta(alpha, beta);
            }
            const double total_count = n + m;
            const double alpha_plus_beta = alpha + beta;
            const double total = alpha_plus_beta + total_count;
            // the removed counts are integers, so their sums are exact
            const double removed_plus = stan_log_beta_ratio_removed(alpha, n)
                                        + stan_log_beta_ratio_removed(beta, m);
            const double removed_minus
                = stan_log_beta_ratio_removed(alpha_plus_beta, total_count);
            const double term_alpha = stan_log_beta_ratio_term(alpha, n);
            const double term_beta = stan_log_beta_ratio_term(beta, m);
            const double term_total
                = stan_log_beta_ratio_term(alpha_plus_beta, total_count);
            // n log(p) + m log(q) with p + q = 1: take the log1m form for
            // the larger of p and q, so that a ratio near 1 keeps its digits
            const double p = (alpha + n) / total;
            const double q = (beta + m) / total;
            const double log_n = (p < q) ? n * log(p) : n * log1p(-q);
            const double log_m = (p < q) ? m * log1p(-p) : m * log(q);
            const double size_stirling = fabs(term_alpha) + fabs(term_beta)
                                         + fabs(term_total) + fabs(log_n)
                                         + fabs(log_m);
            const double size_lbeta
                = stan_lbeta_terms_size(alpha, beta)
                  + stan_lbeta_terms_size(alpha + n, beta + m);
            if (size_stirling <= size_lbeta) {
              const double denominator
                  = lgamma_stirling_diff(alpha) + lgamma_stirling_diff(beta)
                    - lgamma_stirling_diff(alpha_plus_beta);
              return term_alpha + term_beta - term_total
                     + (removed_plus - removed_minus) + log_n + log_m
                     + lgamma_stirling_diff(alpha + n)
                     + lgamma_stirling_diff(beta + m)
                     - lgamma_stirling_diff(total) - denominator;
            }
            return stan_lbeta(alpha + n, beta + m) - stan_lbeta(alpha, beta);
          }
          // \cond
          ) "\n#endif\n";  // NOLINT
// \endcond

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
