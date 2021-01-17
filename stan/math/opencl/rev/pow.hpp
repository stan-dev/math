#ifndef STAN_MATH_OPENCL_REV_POW_HPP
#define STAN_MATH_OPENCL_REV_POW_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the first argument raised to the power of the second
 * argument.
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
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = pow(value_of(a_arena), value_of(b_arena));

  return make_callback_var(
      res_val,
      [a_arena, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto isnan_res = isnan(res.val());
        auto constant_nan = constant(NOT_A_NUMBER, res.rows(), res.cols());
        auto zero_a_arena = value_of(a_arena) == 0.0;
        auto zeros = constant(0, res.rows(), res.cols());
        auto a_deriv = select(
            isnan_res, constant_nan,
            select(zero_a_arena, zeros,
                   elt_multiply(res.adj(),
                                elt_multiply(value_of(b_arena),
                                             elt_divide(res.val(),
                                                        value_of(a_arena))))));
        auto b_deriv = select(
            isnan_res, constant_nan,
            select(zero_a_arena, zeros,
                   elt_multiply(res.adj(), elt_multiply(log(value_of(a_arena)),
                                                        res.val()))));

        results(adjoint_of(a_arena), adjoint_of(b_arena))
            += expressions(calc_if<is_var<T_a>::value>(a_deriv),
                           calc_if<is_var<T_b>::value>(b_deriv));
      });
}

/**
 * Return the first argument raised to the power of the second
 * argument.
 *
 * @tparam T_a type of first expression of the base
 * @tparam T_b type of scalar exponent
 * @param a first expression of base
 * @param b scalar exponent
 * @return base raised to the power of the exponent
 */
template <typename T_a, typename T_b,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_a>* = nullptr,
          require_stan_scalar_t<T_b>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> pow(T_a&& a, const T_b& b) {
  const arena_t<T_a>& a_arena = std::forward<T_a>(a);

  matrix_cl<double> res_val = pow(value_of(a_arena), value_of(b));

  return make_callback_var(
      res_val, [a_arena, b](const vari_value<matrix_cl<double>>& res) mutable {
        auto isnan_res = isnan(res.val());
        auto constant_nan = constant(NOT_A_NUMBER, res.rows(), res.cols());
        auto zero_a_arena = value_of(a_arena) == 0.0;
        auto zeros = constant(0, res.rows(), res.cols());
        if (!is_constant<T_a>::value) {
          auto& a_adj = forward_as<var_value<matrix_cl<double>>>(a_arena).adj();
          a_adj
              = a_adj
                + select(
                      isnan_res, constant_nan,
                      select(zero_a_arena, zeros,
                             elt_multiply(
                                 res.adj(),
                                 elt_multiply(value_of(b),
                                              elt_divide(res.val(),
                                                         value_of(a_arena))))));
        }
        if (!is_constant<T_b>::value) {
          adjoint_of(b) += sum(select(
              isnan_res, constant_nan,
              select(
                  zero_a_arena, zeros,
                  elt_multiply(res.adj(), elt_multiply(log(value_of(a_arena)),
                                                       res.val())))));
        }
      });
}

/**
 * Return the first argument raised to the power of the second
 * argument.
 *
 * @tparam T_a type of scalar exponent
 * @tparam T_b type of first expression of the base
 * @param a scalar base
 * @param b first expression of exponent
 * @return base raised to the power of the exponent
 */
template <typename T_a, typename T_b,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_b>* = nullptr,
          require_stan_scalar_t<T_a>* = nullptr,
          require_any_var_t<T_a, T_b>* = nullptr>
inline var_value<matrix_cl<double>> pow(const T_a& a, T_b&& b) {
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  matrix_cl<double> res_val = pow(value_of(a), value_of(b_arena));

  return make_callback_var(
      res_val, [a, b_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto isnan_res = isnan(res.val());
        auto constant_nan = constant(NOT_A_NUMBER, res.rows(), res.cols());
        auto zero_a_arena = value_of(a) == 0.0;
        auto zeros = constant(0, res.rows(), res.cols());
        if (!is_constant<T_a>::value) {
          adjoint_of(a) += sum(select(
              isnan_res, constant_nan,
              select(zero_a_arena, zeros,
                     elt_multiply(
                         res.adj(),
                         elt_multiply(value_of(b_arena),
                                      elt_divide(res.val(), value_of(a)))))));
        }
        if (!is_constant<T_b>::value) {
          auto& b_adj = forward_as<var_value<matrix_cl<double>>>(b_arena).adj();
          b_adj = b_adj
                  + select(isnan_res, constant_nan,
                           select(zero_a_arena, zeros,
                                  elt_multiply(res.adj(),
                                               elt_multiply(log(value_of(a)),
                                                            res.val()))));
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
