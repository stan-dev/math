#ifndef STAN_MATH_OPENCL_PRIM_LOG_SUM_EXP_HPP
#define STAN_MATH_OPENCL_PRIM_LOG_SUM_EXP_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
#include <stan/math/prim/fun/constants.hpp>

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
 * @return the log of the sum of the exponentiated vector values
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline double log_sum_exp(const T& a) {
  using std::log;
  if (a.size() == 0) {
    return NEGATIVE_INFTY;
  }
  if (stan::internal::is_trivial_kg_expression<T>::value) {
    double a_max = from_matrix_cl(max_2d(a)).maxCoeff();
    if (!std::isfinite(a_max)) {
      return a_max;
    }
    return a_max + log(sum(exp(a - a_max)));
  } else {
    matrix_cl<double> a_eval;
    matrix_cl<double> a_max_cl;
    results(a_eval, a_max_cl) = expressions(a, max_2d(a));
    double a_max = from_matrix_cl(a_max_cl).maxCoeff();
    if (!std::isfinite(a_max)) {
      return a_max;
    }
    return a_max + log(sum(exp(a_eval - a_max)));
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
