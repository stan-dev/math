#ifndef STAN_MATH_REV_FUN_BETA_HPP
#define STAN_MATH_REV_FUN_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/grad_2F1.hpp>

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
  using inner_return_t = decltype(beta(value_of(a), value_of(b)));
  using return_t = return_var_matrix_t<inner_return_t, T1, T2>;
  arena_t<ref_type_t<T1>> arena_a = a;
  arena_t<ref_type_t<T2>> arena_b = b;

  return_t res = beta(value_of(arena_a), value_of(arena_b));
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
  return res;
}

template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr,
          require_any_var_t<T1, T2, T3>* = nullptr>
inline var beta(const T1 a, const T2 b, const T3 x) {
  double a_val = value_of(a);
  double b_val = value_of(b);
  double x_val = value_of(x);
  double res_val = beta(a_val, b_val, x_val);
  return make_callback_var(
    res_val, [a, b, x, a_val, b_val, x_val](const auto& vi) mutable {
      double log_x = log(x_val);
      double log1m_x = log1m(x_val);

      if constexpr (!is_constant_all<T1, T2>::value) {
        auto grad_tuple = grad_2F1(1.0, a + b, 1.0 + a, x_val);
        double grad_mult = exp(b_val * log1m_x + a_val * log_x - log(a_val));

        if constexpr (!is_constant<T1>::value) {
          double da_grads = std::get<1>(grad_tuple) + std::get<2>(grad_tuple);
          double a_adj = vi.val() * (log_x - inv(a_val)) + da_grads * grad_mult;
          a.adj() += vi.adj() * a_adj;
        }

        if constexpr (!is_constant<T2>::value) {
          double db_grad = std::get<1>(grad_tuple);
          double b_adj = vi.val() * log1m_x + db_grad * grad_mult;
          b.adj() += vi.adj() * b_adj;
        }
      }

      if constexpr (!is_constant<T3>::value) {
        double x_adj = exp((b_val - 1.0) * log1m_x + (a_val - 1.0) * log_x);
        x.adj() += vi.adj() * x_adj;
      }
    });
}

}  // namespace math
}  // namespace stan
#endif
