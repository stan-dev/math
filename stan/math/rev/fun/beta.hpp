#ifndef STAN_MATH_REV_FUN_BETA_HPP
#define STAN_MATH_REV_FUN_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/prim/fun/beta.hpp>

namespace stan {
namespace math {

/**
 * Returns the beta function and gradients for two var inputs.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \left(\psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]

   \f[
    \frac{\partial }{\partial b} = \left(\psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @param a var Argument
 * @param b var Argument
 * @return Result of beta function
 */
template <typename T1, typename T2,
          require_all_not_std_vector_t<T1, T2>* = nullptr,
          require_return_type_t<is_var, T1, T2>* = nullptr>
inline auto beta(const T1& a, const T2& b) {
  arena_t<ref_type_t<T1>> arena_a = a;
  arena_t<ref_type_t<T2>> arena_b = b;

  const auto& beta_val = beta(value_of(arena_a), value_of(arena_b));
  using return_type_t = return_var_matrix_t<decltype(beta_val), T1, T2>;
  arena_t<return_type_t> res(beta_val);

  reverse_pass_callback([arena_a, arena_b, res]() mutable {
    auto&& a_array = as_array_or_scalar(arena_a);
    auto&& b_array = as_array_or_scalar(arena_b);
    const auto& res_array = as_array_or_scalar(res);
    const auto& digamma_ab = digamma(value_of(a_array) + value_of(b_array));
    const auto& adj_val = res_array.adj() * res_array.val();

    if constexpr (!is_constant<T1>::value) {
      const auto& a_adj = adj_val * (digamma(a_array.val()) - digamma_ab);
      if constexpr (is_stan_scalar<T1>::value) {
        a_array.adj() += sum(a_adj);
      } else {
        a_array.adj() += a_adj;
      }
    }
    if constexpr (!is_constant<T2>::value) {
      const auto& b_adj = adj_val * (digamma(b_array.val()) - digamma_ab);
      if constexpr (is_stan_scalar<T2>::value) {
        b_array.adj() += sum(b_adj);
      } else {
        b_array.adj() += b_adj;
      }
    }
  });
  return return_type_t(res);
}

}  // namespace math
}  // namespace stan
#endif
