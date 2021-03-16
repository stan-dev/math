#ifndef STAN_MATH_OPENCL_REV_SOFTMAX_HPP
#define STAN_MATH_OPENCL_REV_SOFTMAX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/log_sum_exp.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns softmax of given argument.
 *
 * @tparam T type of the argument
 *
 * @param A argument
 * @return Softmax of the argument
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> softmax(const var_value<T>& A) {
  if (A.size() == 0) {
    return A;
  }
  return make_callback_var(
      softmax(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += elt_multiply(
            res.val(), (res.adj() - dot_product(res.adj(), res.val())));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
