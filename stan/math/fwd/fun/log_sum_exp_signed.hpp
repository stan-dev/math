#ifndef STAN_MATH_FWD_FUN_LOG_SUM_EXP_SIGNED_HPP
#define STAN_MATH_FWD_FUN_LOG_SUM_EXP_SIGNED_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log_sum_exp_signed.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.  The matrix may be a full matrix, a vector,
 * or a row vector. Additionally, a matching container of 'signs' indicates
 * whether the exponentiated input should be added or substracted.
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
template <typename T1, typename T2,
          require_container_st<is_fvar, T1>* = nullptr,
          require_container_st<std::is_integral, T2>* = nullptr>
inline auto log_sum_exp_signed(const T1& x, const T2& signs) {
  return apply_vector_unary<ref_type_t<decltype(to_vector(x))>>::reduce(
      to_ref(to_vector(x)), [&](const auto& v) {
        using T_fvar_inner = typename value_type_t<decltype(v)>::Scalar;
        using vec_type = Eigen::Matrix<T_fvar_inner, -1, 1>;
        Eigen::Map<const Eigen::VectorXi> int_vec_map(signs.data(),
                                                      signs.size());
        vec_type vals = v.val();
        vec_type exp_vals = vals.array().exp().matrix()
                                              .cwiseProduct(int_vec_map);

        return fvar<T_fvar_inner>(
            log_sum_exp_signed(vals, int_vec_map),
            v.d().cwiseProduct(exp_vals).sum() / exp_vals.sum());
      });
}

}  // namespace math
}  // namespace stan
#endif
