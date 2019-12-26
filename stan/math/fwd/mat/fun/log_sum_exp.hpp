#ifndef STAN_MATH_FWD_MAT_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_FWD_MAT_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_vector_unary.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/mat/fun/log_sum_exp.hpp>

namespace stan {
namespace math {

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.  The matrix may be a full matrix, a vector,
 * a row vector, or a container of these.
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @tparam T Type of input vector or matrix.
 * @param[in] x Matrix of specified values.
 * @return The log of the sum of the exponentiated vector values.
 */
template <typename T, require_t<is_fvar<scalar_type_t<T>>>...>
inline auto log_sum_exp(T&& x) {
  return apply_vector_unary<T>::reduce(std::forward<T>(x), [&](auto& v){
    using T_fvar_inner = typename value_type_t<decltype(v)>::Scalar;
    using mat_type = Eigen::Matrix<T_fvar_inner, -1, -1>;
    mat_type vals = v.val();
    mat_type exp_vals = vals.array().exp();

    return fvar<T_fvar_inner>(log_sum_exp(vals),
                   v.d().cwiseProduct(exp_vals).sum() / exp_vals.sum());
  });
}

}  // namespace math
}  // namespace stan
#endif
