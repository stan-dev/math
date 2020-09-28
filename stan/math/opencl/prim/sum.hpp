#ifndef STAN_MATH_OPENCL_PRIM_SUM_HPP
#define STAN_MATH_OPENCL_PRIM_SUM_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Calculates sum of given kernel generator expression.
 * @tparam T type of the expression
 * @param m expression to sum
 * @return sum of given expression
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
value_type_t<T> sum(const T& m) {
  matrix_cl<value_type_t<T>> res;
  if (m.rows() == 1) {
    // using rowwise_sum would run just 1 thread
    res = colwise_sum(transpose(m));
  } else if (m.cols() == 1) {
    res = colwise_sum(m);
  } else {
    res = colwise_sum(rowwise_sum(m));
  }
  return sum(from_matrix_cl<Eigen::Dynamic, 1>(res));
}

}  // namespace math
}  // namespace stan

#endif
#endif
