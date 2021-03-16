#ifndef STAN_MATH_OPENCL_PRIM_SOFTMAX_HPP
#define STAN_MATH_OPENCL_PRIM_SOFTMAX_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified vector.
 *
 * @tparam T type of the vector
 * @param[in] a Vector to transform.
 * @return softmax of the vector
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline matrix_cl<double> softmax(const T& a) {
  check_vector("softmax (OpenCL)", "a", a);
  if (a.size() == 0) {
    return a;
  }
  matrix_cl<double> theta;
  if (stan::internal::is_trivial_kg_expression<T>::value) {
    matrix_cl<double> a_max = max_2d(a);
    theta = exp(a - from_matrix_cl(a_max).maxCoeff());
  } else {
    matrix_cl<double> a_eval;
    matrix_cl<double> a_max;
    results(a_eval, a_max) = expressions(a, max_2d(a));
    theta = exp(a_eval - from_matrix_cl(a_max).maxCoeff());
  }
  return elt_divide(theta, sum(theta));
}

}  // namespace math
}  // namespace stan

#endif
#endif
