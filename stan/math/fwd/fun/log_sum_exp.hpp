#ifndef STAN_MATH_FWD_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_FWD_FUN_LOG_SUM_EXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> log_sum_exp(const fvar<T>& x1, const fvar<T>& x2) {
  using std::exp;
  return fvar<T>(log_sum_exp(x1.val_, x2.val_),
                 x1.d_ / (1 + exp(x2.val_ - x1.val_))
                     + x2.d_ / (exp(x1.val_ - x2.val_) + 1));
}

template <typename T>
inline fvar<T> log_sum_exp(double x1, const fvar<T>& x2) {
  using std::exp;
  if (x1 == NEGATIVE_INFTY) {
    return fvar<T>(x2.val_, x2.d_);
  }
  return fvar<T>(log_sum_exp(x1, x2.val_), x2.d_ / (exp(x1 - x2.val_) + 1));
}

template <typename T>
inline fvar<T> log_sum_exp(const fvar<T>& x1, double x2) {
  return log_sum_exp(x2, x1);
}

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
template <typename T, require_container_st<is_fvar, T>* = nullptr>
inline auto log_sum_exp(const T& x) {
  return apply_vector_unary<ref_type_t<T>>::reduce(
      to_ref(x), [&](const auto& v) {
        using T_fvar_inner = typename value_type_t<decltype(v)>::Scalar;
        using mat_type = Eigen::Matrix<T_fvar_inner, -1, -1>;
        mat_type vals = v.val();
        mat_type exp_vals = vals.array().exp();

        return fvar<T_fvar_inner>(
            log_sum_exp(vals),
            v.d().cwiseProduct(exp_vals).sum() / exp_vals.sum());
      });
}

}  // namespace math
}  // namespace stan
#endif
