#ifndef STAN_MATH_OPENCL_REV_SQUARED_DISTANCE_HPP
#define STAN_MATH_OPENCL_REV_SQUARED_DISTANCE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/prim/squared_distance.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Returns the squared distance.
 *
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return the squared distance of the input.
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<double> squared_distance(T_a&& a, T_b&& b) {
  arena_t<T_a> a_arena = std::forward<T_a>(a);
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  return make_callback_var(
      squared_distance(value_of(a_arena), value_of(b_arena)),
      [a_arena, b_arena](const vari_value<double>& res) mutable {
        auto res_two_mult_diff = elt_multiply(
            res.adj(), 2.0 * (value_of(a_arena) - value_of(b_arena)));
        adjoint_results(a_arena, b_arena)
            += expressions(res_two_mult_diff, -res_two_mult_diff);
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
