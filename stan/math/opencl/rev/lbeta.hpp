#ifndef STAN_MATH_OPENCL_REV_LBETA_HPP
#define STAN_MATH_OPENCL_REV_LBETA_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise `lbeta()` on two input kernel
 * generator expression
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return elementwise `lbeta()`
 */
template <typename T_a, typename T_b,
          require_all_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr,
          require_any_not_stan_scalar_t<T_a, T_b>* = nullptr>
inline auto lbeta(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = std::forward<T_a>(a);
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  return make_callback_var(
      lbeta(value_of(a_arena), value_of(b_arena)),
      [a_arena, b_arena](vari_value<matrix_cl<double>>& res) mutable {
        auto digamma_ab = digamma(value_of(a_arena) + value_of(b_arena));
        adjoint_results(a_arena, b_arena) += expressions(
            elt_multiply(res.adj(), digamma(value_of(a_arena)) - digamma_ab),
            elt_multiply(res.adj(), digamma(value_of(b_arena)) - digamma_ab));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
