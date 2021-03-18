#ifndef STAN_MATH_OPENCL_REV_LOG_SOFTMAX_HPP
#define STAN_MATH_OPENCL_REV_LOG_SOFTMAX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/log_softmax.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return log of the softmax of the specified vector.
 *
 * @tparam T type of the vector
 *
 * @param[in] A vector
 * @return log of the softmax
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> log_softmax(const var_value<T>& A) {
  return make_callback_var(
      log_softmax(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += res.adj() - sum(res.adj()) * exp(res.val());
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
