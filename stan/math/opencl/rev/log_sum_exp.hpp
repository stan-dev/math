#ifndef STAN_MATH_OPENCL_REV_LOG_SUM_EXP_HPP
#define STAN_MATH_OPENCL_REV_LOG_SUM_EXP_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/log_sum_exp.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.
 *
 * @tparam T type of the vector
 *
 * @param[in] A vector
 * @return the log of the sum of the exponentiated vector values
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var log_sum_exp(const var_value<T>& A) {
  return make_callback_var(log_sum_exp(A.val()), [A](vari& res) mutable {
    A.adj() += res.adj() * exp(A.val() - res.val());
  });
}

}  // namespace math
}  // namespace stan

#endif
#endif
