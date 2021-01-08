#ifndef STAN_MATH_OPENCL_REV_POW_HPP
#define STAN_MATH_OPENCL_REV_POW_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the first argument raised to the power of the second
 * argument.  At least one of the arguments must be a complex number.
 *
 * @tparam T_a type of first expression of the base
 * @tparam T_b type of second expression of the exponent
 * @param a first expression of base
 * @param b second expression of exponent
 * @return base raised to the power of the exponent
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> pow(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = std::forward<T_a>(a);
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = pow(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto a_deriv = select(
            isnan(value_of(res)),
            constant(NOT_A_NUMBER, res.rows(), res.cols()),
            select(
                value_of(a_arena) == 0.0, constant(0, res.rows(), res.cols()),
                elt_multiply(
                    res.adj(),
                    elt_multiply(value_of(b_arena),
                                 elt_divide(res.val(), value_of(a_arena))))));
        auto b_deriv = select(
            isnan(value_of(res)),
            constant(NOT_A_NUMBER, res.rows(), res.cols()),
            select(
                value_of(a_arena) == 0.0, constant(0, res.rows(), res.cols()),
                elt_multiply(res.adj(),
                             elt_multiply(log(value_of(a_arena)), res.val()))));

        results(adjoint_of(a_arena), adjoint_of(b_arena))
            += expressions(calc_if<is_var<T_a>::value>(a_deriv),
                           calc_if<is_var<T_b>::value>(b_deriv));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
