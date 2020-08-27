#ifndef STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#define STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2,
          require_all_prim_or_rev_kernel_expressions<T1, T2>* = nullptr>
auto multiply(const T1& a, const T2& b) {
  check_size_match("multiply ((OpenCL))", "A.cols()", A.cols(), "B.rows()",
                   B.rows());
  const auto& a_eval = eval(a);  // TODO check eval, arena types
  const auto& b_eval = eval(b);

  var_value<matrix_cl<double>> res;
  res.val() = value_of(a_eval) * value_of(b_eval);

  reverse_pass_callback([a_eval, b_eval, res]() mutable {
    if (!is_constant<T_a>::value) {
      a_eval.adj() = a_eval.adj() + res.adj() * transpose(value_of(b_eval));
    }
    if (!is_constant<T_b>::value) {
      b_eval.adj() = b_eval.adj() + transpose(value_of(a_eval)) * res.adj();
    }
  });
  return res;
}

}  // namespace math

#endif
#endif
