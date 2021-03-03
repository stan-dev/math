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
 * Returns the dot product.
 *
 * @tparam T1 type of the first vector
 * @tparam T2 type of the second vector
 *
 * @param[in] v1 First vector.
 * @param[in] v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error if sizes of v1 and v2 do not match.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> log_softmax(const var_value<T>& A) {
  return make_callback_var(
      log_softmax(A.val()),
      [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += res.adj() - sum(res.adj()) * exp(res.val());
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
