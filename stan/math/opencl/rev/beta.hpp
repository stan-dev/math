#ifndef STAN_MATH_OPENCL_REV_BETA_HPP
#define STAN_MATH_OPENCL_REV_BETA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise `beta()` on two input kernel
 * generator expression
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return elementwise `beta()`
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline auto beta(T_a&& a, T_b&& b) {
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  var_value<matrix_cl<double>> res = beta(value_of(a_arena), value_of(b_arena));

  reverse_pass_callback([a_arena, b_arena, res]() mutable {
    auto adj_val = elt_multiply(res.adj(), res.val());
    auto digamma_ab = digamma(value_of(a_arena) + value_of(b_arena));
    adjoint_results(a_arena, b_arena) += expressions(
        elt_multiply(adj_val, (digamma(value_of(a_arena)) - digamma_ab)),
        elt_multiply(adj_val, (digamma(value_of(b_arena)) - digamma_ab)));
  });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
