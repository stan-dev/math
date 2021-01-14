#ifndef STAN_MATH_OPENCL_REV_LDEXP_HPP
#define STAN_MATH_OPENCL_REV_LDEXP_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the elementwise `ldexp()` of the input
 * `var_value<matrix_cl<double>>` and kernel generator expression.
 *
 * @param a input rev kernel generator expression representing
 * significands
 * @param b input kernel generator expression representing
 * the integer exponents.
 * @return Elementwise `ldexp()` of the input argument.
 */
template <typename T_a, typename T_b,
          require_all_kernel_expressions_and_none_scalar_t<T_a>* = nullptr,
          require_all_kernel_expressions_t<T_b>* = nullptr>
inline var_value<matrix_cl<double>> ldexp(const var_value<T_a>& a, T_b&& b) {
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  var_value<matrix_cl<double>> res = ldexp(a.val(), b);

  reverse_pass_callback([a, b_arena, res]() mutable {
    a.adj() = a.adj() + ldexp(res.adj(), value_of(b_arena));
  });

  return res;
}

/**
 * Returns the elementwise `ldexp()` of the input
 * `var_value<double>` and kernel generator expression.
 *
 * @param a input rev kernel generator expression representing
 * significands
 * @param b input kernel generator expression representing
 * the integer exponents.
 * @return Elementwise `ldexp()` of the input argument.
 */
template <typename T_b,
          require_all_kernel_expressions_and_none_scalar_t<T_b>* = nullptr>
inline var_value<matrix_cl<double>> ldexp(const var_value<double>& a, T_b&& b) {
  const arena_t<T_b>& b_arena = std::forward<T_b>(b);

  var_value<matrix_cl<double>> res = ldexp(a.val(), b);

  reverse_pass_callback([a, b_arena, res]() mutable {
    a.adj() = a.adj() + sum(ldexp(res.adj(), value_of(b_arena)));
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
