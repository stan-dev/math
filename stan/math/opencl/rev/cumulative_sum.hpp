#ifndef STAN_MATH_OPENCL_REV_CUMULATIVE_SUM_HPP
#define STAN_MATH_OPENCL_REV_CUMULATIVE_SUM_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/cumulative_sum.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the cumulative sum of the specified vector.
 *
 * The cumulative sum of a vector of values \code{x} is the
 *
 * \code x[0], x[1] + x[2], ..., x[1] + , ..., + x[x.size()-1] @endcode
 *
 * @tparam T scalar type of the vector
 * @param A Vector of values
 * @return Cumulative sum of values
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> cumulative_sum(const var_value<T>& A) {
  return make_callback_var(
      cumulative_sum(A.val()), [A](vari_value<matrix_cl<double>>& res) mutable {
        A.adj() += reverse(cumulative_sum(reverse(res.adj())));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
