#ifndef STAN_MATH_PRIM_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Calculates the log sum of exponentials without overflow.
 *
 * \f$\log (\exp(a) + \exp(b)) = m + \log(\exp(a-m) + \exp(b-m))\f$,
 *
 * where \f$m = max(a, b)\f$.
 *
   \f[
   \mbox{log\_sum\_exp}(x, y) =
   \begin{cases}
     \ln(\exp(x)+\exp(y)) & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial x} =
   \begin{cases}
     \frac{\exp(x)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial y} =
   \begin{cases}
     \frac{\exp(y)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T1 type of the first variable
 * @tparam T2 type of the second variable
 * @param a the first variable
 * @param b the second variable
 */
template <typename T1, typename T2, require_all_not_st_var<T1, T2>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_sum_exp(const T2& a, const T1& b) {
  if (a == NEGATIVE_INFTY) {
    return b;
  }
  if (a == INFTY && b == INFTY) {
    return INFTY;
  }
  if (a > b) {
    return a + log1p_exp(b - a);
  }
  return b + log1p_exp(a - b);
}

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.  The matrix may be a full matrix, a vector,
 * a row vector, or a container of these.
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @tparam T type of input vector or matrix
 * @param[in] x matrix of specified values
 * @return The log of the sum of the exponentiated vector values.
 */
template <typename T, require_container_st<std::is_arithmetic, T>* = nullptr>
inline auto log_sum_exp(const T& x) {
  return apply_vector_unary<T>::reduce(x, [&](const auto& v) {
    if (v.size() == 0) {
      return NEGATIVE_INFTY;
    }
    const auto& v_ref = to_ref(v);
    const double max = v_ref.maxCoeff();
    if (!std::isfinite(max)) {
      return max;
    }
    return max + std::log((v_ref.array() - max).exp().sum());
  });
}

/**
 * Enables the vectorized application of the log_sum_exp function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return log_sum_exp function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_sum_exp(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [](const auto& c, const auto& d) { return log_sum_exp(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
