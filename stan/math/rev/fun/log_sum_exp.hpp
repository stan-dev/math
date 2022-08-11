#ifndef STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the log of the sum of the exponentiation of the arguments.
 *
 * @param[in] a first argument
 * @param[in] b second argument
 * @return log-sum-exp of arguments
 */
inline var log_sum_exp(const var& a, const var& b) {
  return make_callback_var(log_sum_exp(a.val(), b.val()),
			   [a, b](const auto& y) mutable {
			     a.adj() += y.adj() * inv_logit(a.val() - b.val());
			     b.adj() += y.adj() * inv_logit(b.val() - a.val());
			   });
}

/**
 * Returns the log of the sum of the exponentiation of the arguments.
 *
 * @param[in] a first argument
 * @param[in] b second argument
 * @return log-sum-exp of arguments
 */
inline var log_sum_exp(const var& a, double b) {
  return make_callback_var(log_sum_exp(a.val(), b),
			   [a, b](const auto& y) mutable {
			     a.adj() += a.val() == NEGATIVE_INFTY
			       ? y.adj()
			       : y.adj() * inv_logit(a.val() - b);
			   });
}

/**
 * Returns the log of the sum of the exponentiation of the arguments.
 *
 * @param[in] a first argument
 * @param[in] b second argument
 * @return log-sum-exp of arguments
 */
inline var log_sum_exp(double a, const var& b) {
  return log_sum_exp(b, a);
}

/**
 * Returns the log of the sum of the input exponentiated.
 *
 * @tparam T container type inheriting from `EigenBase` with scalar type `var`
 * @param v argument container
 * @return log-sum-exp of container elements
 */
template <typename T, require_eigen_st<is_var, T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline var log_sum_exp(const T& v) {
  arena_t<decltype(v)> arena_v = v;
  return make_callback_var(log_sum_exp(v.val()),
			   [arena_v](const auto& res) mutable {
			     arena_v.adj()
			       += res.adj() * (arena_v.val().array().val() - res.val()).exp().matrix();
  });		   
}

/**
 * Returns the log of the sum of the exponentiated argument.
 *
 * @tparam T matrix or vector type of argument
 * @param x argument
 * @return log-sum-exp of argument
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var log_sum_exp(const T& x) {
  return make_callback_var(log_sum_exp(x.val()), [x](const auto& res) mutable {
    x.adj() += res.adj() * (x.val().array().val() - res.val()).exp().matrix();
  });
}

/**
 * Returns the log of the sum of the exponentiated argument.
 *
 * @tparam T type of container (standard vector)
 * @param x argument container
 * @return log-sum-exp of elements of argument
 */
template <typename T, require_std_vector_st<is_var, T>* = nullptr>
inline auto log_sum_exp(const T& x) {
  return apply_vector_unary<T>::reduce(
      x, [](const auto& v) { return log_sum_exp(v); });
}

}  // namespace math
}  // namespace stan
#endif
