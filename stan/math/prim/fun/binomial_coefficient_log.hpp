#ifndef STAN_MATH_PRIM_FUN_BINOMIAL_COEFFICIENT_LOG_HPP
#define STAN_MATH_PRIM_FUN_BINOMIAL_COEFFICIENT_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

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
 *
   \f[
   \mbox{binomial\_coefficient\_log}(x, y) =
   \begin{cases}
     \textrm{error} & \mbox{if } y > x + 1 \textrm{ or } y < -1 \textrm{ or } x
 < -1\\
     \ln\Gamma(x+1) & \mbox{if } -1 < y < x + 1 \\
     \quad -\ln\Gamma(y+1)& \\
     \quad -\ln\Gamma(x-y+1)& \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{binomial\_coefficient\_log}(x, y)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } y > x + 1 \textrm{ or } y < -1 \textrm{ or } x
 < -1\\
     \Psi(x+1) & \mbox{if } 0\leq y \leq x \\
     \quad -\Psi(x-y+1)& \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{binomial\_coefficient\_log}(x, y)}{\partial y} =
   \begin{cases}
     \textrm{error} & \mbox{if } y > x + 1 \textrm{ or } y < -1 \textrm{ or } x
 < -1\\
     -\Psi(y+1) & \mbox{if } 0\leq y \leq x \\
     \quad +\Psi(x-y+1)& \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 *  This function is numerically more stable than naive evaluation via lgamma.
 *
 * @tparam T_N type of N.
 * @tparam T_n type of n.
 * @param N total number of objects.
 * @param n number of objects chosen.
 * @return log (N choose n).
 */

template <typename T_N, typename T_n>
inline return_type_t<T_N, T_n> binomial_coefficient_log(const T_N N,
                                                        const T_n n) {
  if (is_any_nan(N, n)) {
    return stan::math::NOT_A_NUMBER;
  }

  const T_N N_plus_1 = N + 1;

  static const char* function = "binomial_coefficient_log";
  check_greater_or_equal(function, "first argument", N, -1);
  check_greater_or_equal(function, "second argument", n, -1);
  check_greater_or_equal(function, "(first argument - second argument + 1)",
                         N_plus_1 - n, 0.0);

  if (n == 0) {
    return 0;
  }
  if (N / 2 < n) {
    return binomial_coefficient_log(N, N - n);
  } 
  if (N_plus_1 < lgamma_stirling_diff_useful) {
    return lgamma(N_plus_1) - lgamma(n + 1) - lgamma(N_plus_1 - n);
  } 
  return -lbeta(N - n + 1, n + 1) - log(N_plus_1);
}

}  // namespace math
}  // namespace stan
#endif
