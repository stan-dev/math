#ifndef STAN_MATH_PRIM_FUN_MULTIPLY_LOG_HPP
#define STAN_MATH_PRIM_FUN_MULTIPLY_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

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
     x\ln y & \mbox{if } x, y\neq 0 \\[6pt]
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{multiply\_log}(x, y)}{\partial x} =
   \begin{cases}
     \ln y \\[6pt]
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{multiply\_log}(x, y)}{\partial y} =
   \begin{cases}
     \frac{x}{y} \\[6pt]
   \end{cases}
   \f]
 *
 * @tparam T_a type of the first variable
 * @tparam T_b type of the second variable
 * @param a the first variable
 * @param b the second variable
 * @return a * log(b)
 */
template <typename T_a, typename T_b,
          require_all_arithmetic_t<T_a, T_b>* = nullptr>
inline return_type_t<T_a, T_b> multiply_log(const T_a a, const T_b b) {
  using std::log;
  if (b == 0.0 && a == 0.0) {
    return 0.0;
  }

  return a * log(b);
}

/**
 * Enables the vectorised application of the multiply_log
 * function, when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return multiply_log function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_all_not_var_matrix_t<T1, T2>* = nullptr>
inline auto multiply_log(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return multiply_log(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
