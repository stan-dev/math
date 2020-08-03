#ifndef STAN_MATH_PRIM_FUN_BINOMIAL_COEFFICIENT_LOG_HPP
#define STAN_MATH_PRIM_FUN_BINOMIAL_COEFFICIENT_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * Return the log of the binomial coefficient for the specified
 * arguments.
 *
 * The binomial coefficient, \f${n \choose k}\f$, read "n choose k", is
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
 * @tparam T_n type of the first argument
 * @tparam T_k type of the second argument
 *
 * @param n total number of objects.
 * @param k number of objects chosen.
 * @return log (n choose k).
 */

template <typename T_n, typename T_k,
          require_all_stan_scalar_t<T_n, T_k>* = nullptr>
inline return_type_t<T_n, T_k> binomial_coefficient_log(const T_n n,
                                                        const T_k k) {
  using T_partials_return = partials_return_t<T_n, T_k>;

  if (is_any_nan(n, k)) {
    return NOT_A_NUMBER;
  }

  // Choosing the more stable of the symmetric branches
  if (n > -1 && k > value_of_rec(n) / 2.0 + 1e-8) {
    return binomial_coefficient_log(n, n - k);
  }

  const T_partials_return n_dbl = value_of(n);
  const T_partials_return k_dbl = value_of(k);
  const T_partials_return n_plus_1 = n_dbl + 1;
  const T_partials_return n_plus_1_mk = n_plus_1 - k_dbl;

  static const char* function = "binomial_coefficient_log";
  check_greater_or_equal(function, "first argument", n, -1);
  check_greater_or_equal(function, "second argument", k, -1);
  check_greater_or_equal(function, "(first argument - second argument + 1)",
                         n_plus_1_mk, 0.0);

  operands_and_partials<T_n, T_k> ops_partials(n, k);

  T_partials_return value;
  if (k_dbl == 0) {
    value = 0;
  } else if (n_plus_1 < lgamma_stirling_diff_useful) {
    value = lgamma(n_plus_1) - lgamma(k_dbl + 1) - lgamma(n_plus_1_mk);
  } else {
    value = -lbeta(n_plus_1_mk, k_dbl + 1) - log1p(n_dbl);
  }

  if (!is_constant_all<T_n, T_k>::value) {
    // Branching on all the edge cases.
    // In direct computation many of those would be NaN
    // But one-sided limits from within the domain exist, all of the below
    // follows from lim x->0 from above digamma(x) == -Inf
    //
    // Note that we have k < n / 2 (see the first branch in this function)
    // se we can ignore the n == k - 1 edge case.
    T_partials_return digamma_n_plus_1_mk = digamma(n_plus_1_mk);

    if (!is_constant_all<T_n>::value) {
      if (n_dbl == -1.0) {
        if (k_dbl == 0) {
          ops_partials.edge1_.partials_[0] = 0;
        } else {
          ops_partials.edge1_.partials_[0] = NEGATIVE_INFTY;
        }
      } else {
        ops_partials.edge1_.partials_[0]
            = (digamma(n_plus_1) - digamma_n_plus_1_mk);
      }
    }
    if (!is_constant_all<T_k>::value) {
      if (k_dbl == 0 && n_dbl == -1.0) {
        ops_partials.edge2_.partials_[0] = NEGATIVE_INFTY;
      } else if (k_dbl == -1) {
        ops_partials.edge2_.partials_[0] = INFTY;
      } else {
        ops_partials.edge2_.partials_[0]
            = (digamma_n_plus_1_mk - digamma(k_dbl + 1));
      }
    }
  }

  return ops_partials.build(value);
}

/**
 * Enables the vectorised application of the binomial coefficient log function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Binomial coefficient log function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto binomial_coefficient_log(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return binomial_coefficient_log(c, d);
  });
}

}  // namespace math
}  // namespace stan
#endif
