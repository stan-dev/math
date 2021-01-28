#ifndef STAN_MATH_OPENCL_REV_LDEXP_HPP
#define STAN_MATH_OPENCL_REV_LDEXP_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
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
          require_all_kernel_expressions_t<T_a, T_b>* = nullptr,
          require_st_integral<T_b>* = nullptr>
inline var_value<matrix_cl<double>> ldexp(const var_value<T_a>& a, T_b&& b) {
  arena_t<T_b> b_arena = std::forward<T_b>(b);

  return make_callback_var(
      ldexp(a.val(), b),
      [a, b_arena](vari_value<matrix_cl<double>>& res) mutable {
        adjoint_results(a) += expressions(ldexp(res.adj(), b_arena));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
