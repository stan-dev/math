#ifndef STAN_MATH_OPENCL_PRIM_VARIANCE_HPP
#define STAN_MATH_OPENCL_PRIM_VARIANCE_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/log_sum_exp.hpp>
#include <stan/math/opencl/prim/mean.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the sample variance of the var_value matrix
 * Raise domain error if size is not greater than zero.
 *
 * @tparam T input matrix type
 * @param a input
 * @return sample variance of the input
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline double variance(const T& a) {
  check_nonzero_size("variance (OpenCL)", "a", a);
  if (a.size() == 1) {
    return 0.0;
  }
  if (stan::internal::is_trivial_kg_expression<T>::value) {
    return sum(square(a - mean(a))) / (a.size() - 1);
  } else {
    matrix_cl<double> a_eval;
    matrix_cl<double> a_sum;
    results(a_eval, a_sum) = expressions(a, sum_2d(a));
    return sum(square(a_eval - from_matrix_cl(a_sum).sum() / a.size()))
           / (a.size() - 1);
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
