#ifndef STAN_MATH_OPENCL_PRIM_LOG_SUM_EXP_HPP
#define STAN_MATH_OPENCL_PRIM_LOG_SUM_EXP_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/opencl/scalar_cl_reduce.hpp>

namespace stan {
namespace math {

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values. The matrix may be a full matrix, a vector,
 * a row vector.
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @tparam T type of input vector or matrix
 * @param[in] a matrix of specified values
 * @return device scalar holding the log of the sum of the exponentiated
 * vector values
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline opencl::ScalarCl<double> log_sum_exp(const T& a) {
  if (a.size() == 0) {
    return opencl::ScalarCl<double>(scalar_<double>(NEGATIVE_INFTY));
  }
  const auto log_sum_exp_impl
      = [](const auto& x, const opencl::ScalarCl<double>& x_max) {
          opencl::ScalarCl<double> sum_exp = sum(exp(x - x_max));
          auto x_max_op = as_operation_cl(x_max);
          // a non finite maximum is the result; exp(x - x_max) is not used then
          return opencl::ScalarCl<double>(
              select(isfinite(x_max_op),
                     x_max_op + log(as_operation_cl(sum_exp)), x_max_op));
        };
  if constexpr (stan::internal::is_trivial_kg_expression<T>::value) {
    return log_sum_exp_impl(a, opencl::internal::max_scalar_cl(a));
  } else {
    matrix_cl<double> a_eval;
    matrix_cl<double> a_max_partial;
    results(a_eval, a_max_partial) = expressions(a, max_2d(a));
    return log_sum_exp_impl(a_eval,
                            opencl::internal::max_scalar_cl(a_max_partial));
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
