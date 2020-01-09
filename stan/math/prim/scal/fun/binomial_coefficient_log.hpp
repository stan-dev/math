#ifndef STAN_MATH_PRIM_SCAL_FUN_BINOMIAL_COEFFICIENT_LOG_HPP
#define STAN_MATH_PRIM_SCAL_FUN_BINOMIAL_COEFFICIENT_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma_stirling.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/scal/fun/lgamma_stirling_diff.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>

namespace stan {
namespace math {

/**
 * Return the log of the binomial coefficient for the specified
 * arguments.
 *
 * The binomial coefficient, \f${N \choose n}\f$, read "N choose n", is
 * defined for \f$0 \leq n \leq N\f$ by
 *
 * \f${N \choose n} = \frac{N!}{n! (N-n)!}\f$.
 *
 * This function uses Gamma functions to define the log
 * and generalize the arguments to continuous N and n.
 *
 * \f$ \log {N \choose n}
 * = \log \ \Gamma(N+1) - \log \Gamma(n+1) - \log \Gamma(N-n+1)\f$.
 *
   \f[
   \mbox{binomial\_coefficient\_log}(x, y) =
   \begin{cases}
     \textrm{error} & \mbox{if } y > x \textrm{ or } y < 0\\
     \ln\Gamma(x+1) & \mbox{if } 0\leq y \leq x \\
     \quad -\ln\Gamma(y+1)& \\
     \quad -\ln\Gamma(x-y+1)& \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{binomial\_coefficient\_log}(x, y)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } y > x \textrm{ or } y < 0\\
     \Psi(x+1) & \mbox{if } 0\leq y \leq x \\
     \quad -\Psi(x-y+1)& \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{binomial\_coefficient\_log}(x, y)}{\partial y} =
   \begin{cases}
     \textrm{error} & \mbox{if } y > x \textrm{ or } y < 0\\
     -\Psi(y+1) & \mbox{if } 0\leq y \leq x \\
     \quad +\Psi(x-y+1)& \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param N total number of objects.
 * @param n number of objects chosen.
 * @return log (N choose n).
 */

namespace internal {
constexpr double lbeta_stirling_bound = 10;
}

inline double lbeta(const double a, const double b) {
  double x;  // x is the smaller of the two
  double y;
  if (a < b) {
    x = a;
    y = b;
  } else {
    x = b;
    y = a;
  }

  if (y < internal::lbeta_stirling_bound) {
    // both small
    return lgamma(x) + lgamma(y) - lgamma(x + y);
  } else if (x < internal::lbeta_stirling_bound) {
    // y large, x small
    double stirling_diff
        = lgamma_stirling_diff(y) - lgamma_stirling_diff(x + y);
    double log_x_y = log(x + y);  // log_sum_exp(log(x), log(y));
    double stirling = (y - 0.5) * log1p(-x / (x + y)) + x * (1 - log_x_y);
    // std::cout << std::setprecision(18) << x << " " << y << std::endl;
    // std::cout << std::setprecision(18) << lgamma_stirling_diff(y) << " " <<
    // lgamma(y) << " " << lgamma_stirling_diff(x + y) << " " << lgamma(x + y)
    // << std::endl;
    return stirling + lgamma(x) + stirling_diff;
  } else {
    // both large
    double stirling_diff = lgamma_stirling_diff(x) + lgamma_stirling_diff(y)
                           - lgamma_stirling_diff(x + y);
    double stirling = (x - 0.5) * log(x / (x + y)) + y * log1p(-x / (x + y))
                      + 0.5 * (log(2 * stan::math::pi()) - log(y));
    return stirling + stirling_diff;
  }
}

template <typename T_N, typename T_n>
inline return_type_t<T_N, T_n> binomial_coefficient_log(const T_N N,
                                                        const T_n n) {
  const double CUTOFF = 1000;
  const T_N N_plus_1 = N + 1;
  if (n == 0) {
    return 0;
  } else if (N - n < CUTOFF) {
    return lgamma(N_plus_1) - lgamma(n + 1) - lgamma(N_plus_1 - n);
  } else {
    return -lbeta(N - n + 1, n + 1) - log(N_plus_1);
    // return_type_t<T_N, T_n> N_minus_n = N - n;
    // const double one_twelfth = inv(12);
    // return multiply_log(n, N_minus_n) + multiply_log((N + 0.5), N /
    // N_minus_n)
    //        + one_twelfth / N - n - one_twelfth / N_minus_n - lgamma(n + 1);
  }
}

}  // namespace math
}  // namespace stan
#endif
