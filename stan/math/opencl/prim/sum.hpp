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
  if (is_matrix_cl<T>::value && m.size() < 1000) {
    // for small matrices running another kernel is not worth it
    return sum(from_matrix_cl(m));
  }
  matrix_cl<value_type_t<T>> res;
  if (m.rows() <= 8) {
    // without transpose we would use just a few threads in a work group
    res = sum_2d(transpose(m));
  } else {
    res = sum_2d(m);
  }
  return sum(from_matrix_cl(res));
}

}  // namespace math
}  // namespace stan

#endif
#endif
