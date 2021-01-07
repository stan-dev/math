#ifndef STAN_MATH_OPENCL_REV_HYPOT_HPP
#define STAN_MATH_OPENCL_REV_HYPOT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `hypot()` of the input.
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * 
 * @return Elementwise `hypot()` of the input.
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> hypot(T_a&& a, T_b&& b) {
    arena_t<T_a> a_arena = std::forward<T_a>(a);
    arena_t<T_b> b_arena = std::forward<T_b>(b);

    matrix_cl<double> res_val = hypot(value_of(a_arena), value_of(b_arena));

    return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto a_deriv
            = elt_multiply(res.adj(), elt_divide(value_of(a_arena), res.val()));
        auto b_deriv
            = elt_multiply(res.adj(), elt_divide(value_of(b_arena), res.val()));
        results(adjoint_of(a_arena), adjoint_of(b_arena))
            += expressions(calc_if<is_var<T_a>::value>(a_deriv),
                           calc_if<is_var<T_b>::value>(b_deriv));
    });
}

}  // namespace math
}  // namespace stan

#endif
#endif
