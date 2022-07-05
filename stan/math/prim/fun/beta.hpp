#ifndef STAN_MATH_PRIM_FUN_BETA_HPP
#define STAN_MATH_PRIM_FUN_BETA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the beta function applied to the specified
 * arguments.
 *
 * The beta function is defined for \f$a > 0\f$ and \f$b > 0\f$ by
 *
 * \f$\mbox{B}(a, b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}\f$.
 *
   \f[
   \mbox{beta}(\alpha, \beta) =
   \begin{cases}
     \int_0^1 u^{\alpha - 1} (1 - u)^{\beta - 1} \, du & \mbox{if } \alpha,
 \beta>0 \\[6pt] \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta =
 \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \alpha} =
   \begin{cases}
     \left(\psi(\alpha)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \beta} =
   \begin{cases}
     \left(\psi(\beta)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T1 type of first value
 * @tparam T2 type of second value
 * @param a First value
 * @param b Second value
 * @return Beta function applied to the two values.
 */
template <typename T1, typename T2,
          require_all_not_std_vector_t<T1, T2>* = nullptr,
          require_t<
            disjunction<
              is_any_var_matrix<T1, T2>,
              conjunction<is_stan_scalar<T1>,is_stan_scalar<T2>>
            >
          >* = nullptr>
inline auto beta(const T1 a, const T2 b) {
  auto val_fun = [&](auto&& x, auto&& y) {
    using std::exp;
    return exp(subtract(add(lgamma(x), lgamma(y)), lgamma(add(x, y))));
  };

  auto shared_args_fun
      = [&](auto&& val, auto&& adj, auto&& x, auto&& y) {
        auto digamma_xy = digamma(add(x, y));
        return std::tuple<decltype(digamma_xy)>(digamma_xy);
      };

  auto grad_fun_a
    = [&](auto&& val, auto&& adj, auto&& shared_args, auto&& x, auto&& y) {
      const auto& digamma_xy = std::get<0>(shared_args);

      return elt_multiply(val,
                          elt_multiply(adj, subtract(digamma(x), digamma_xy)));
    };
  auto grad_fun_b
    = [&](auto&& val, auto&& adj, auto&& shared_args, auto&& x, auto&& y) {
      const auto& digamma_xy = std::get<0>(shared_args);

      return elt_multiply(val,
                          elt_multiply(adj, subtract(digamma(y), digamma_xy)));
    };
  return function_gradients<true>(
    std::forward_as_tuple(a, b),
    std::forward<decltype(val_fun)>(val_fun),
    std::forward<decltype(shared_args_fun)>(shared_args_fun),
    std::forward_as_tuple(grad_fun_a, grad_fun_b));
}

/**
 * Enables the vectorised application of the beta function, when
 * the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Beta function applied to the two inputs.
 */
template <typename T1, typename T2,
          require_any_container_t<T1, T2>* = nullptr,
          require_all_not_var_matrix_t<T1, T2>* = nullptr>
inline auto beta(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [](const auto& c, const auto& d) { return beta(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
