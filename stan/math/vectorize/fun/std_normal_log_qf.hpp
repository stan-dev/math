
#ifndef STAN_MATH_VECTORIZED_FUN_STD_NORMAL_LOG_QF_HPP
#define STAN_MATH_VECTORIZED_FUN_STD_NORMAL_LOG_QF_HPP
#include <stan/math/prim/prob/std_normal_log_qf.hpp>
namespace stan {
namespace math {

/**
 * A vectorized version of std_normal_log_qf() that accepts std::vectors, Eigen
 * Matrix/Array objects, or expressions, and containers of these.
 *
 * @tparam T type of container
 * @param x container
 * @return inverse unit normal CDF of each value in x
 * @throw std::domain_error if x is not less than or equal to 0
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto std_normal_log_qf(const T& x) {
  return apply_scalar_unary<std_normal_log_qf_fun, T>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

