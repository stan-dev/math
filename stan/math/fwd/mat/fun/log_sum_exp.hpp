#ifndef STAN_MATH_FWD_MAT_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_FWD_MAT_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/mat/fun/log_sum_exp.hpp>

namespace stan {
namespace math {

template <typename T, require_t<is_fvar<scalar_type_t<T>>>...>
inline auto log_sum_exp(T&& x) {
  return apply_vector_unary<T>::reduce(std::forward<T>(x), [](auto& v){
    using T_scalar = value_type_t<T>;
    using T_fvar = typename T_scalar::Scalar;
    using mat_type = Eigen::Matrix<T_fvar, Eigen::Dynamic, Eigen::Dynamic>;
    mat_type vals = v.val();
    mat_type exp_vals = vals.array().exp();

    return fvar<T_fvar>(log_sum_exp(vals),
                   v.d().cwiseProduct(exp_vals).sum() / exp_vals.sum());
  });
}

}  // namespace math
}  // namespace stan
#endif
