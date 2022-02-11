#ifndef STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * @tparam T Scalar.
 * @tparam L Scalar.
 * @param[in] x Unconstrained input
 * @param[in] lb Lower bound
 * @return Constrained matrix
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb) {
  if (unlikely(value_of_rec(lb) == NEGATIVE_INFTY)) {
    return identity_constrain(x, lb);
  } else {
    return add(exp(x), lb);
  }
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * @tparam T Scalar.
 * @tparam L Scalar.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  if (value_of_rec(lb) == NEGATIVE_INFTY) {
    return identity_constrain(x, lb);
  } else {
    lp += x;
    return add(exp(x), lb);
  }
}

/**
 * Specialization of `lb_constrain` to apply a scalar lower bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam L Scalar.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_eigen_t<T>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(T&& x, L&& lb) {
  return eval(x.unaryExpr([lb](auto&& x) { return lb_constrain(x, lb); }));
}

/**
 * Specialization of `lb_constrain` to apply a scalar lower bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam L Scalar.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_eigen_t<T>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  return eval(
      x.unaryExpr([lb, &lp](auto&& xx) { return lb_constrain(xx, lb, lp); }));
}

/**
 * Specialization of `lb_constrain` to apply a matrix of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam L A type inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_eigen_t<T, L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(T&& x, L&& lb) {
  check_matching_dims("lb_constrain", "x", x, "lb", lb);
  return eval(x.binaryExpr(
      lb, [](auto&& x, auto&& lb) { return lb_constrain(x, lb); }));
}

/**
 * Specialization of `lb_constrain` to apply a matrix of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam L A type inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_eigen_t<T, L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  check_matching_dims("lb_constrain", "x", x, "lb", lb);
  return eval(x.binaryExpr(
      lb, [&lp](auto&& xx, auto&& lbb) { return lb_constrain(xx, lbb, lp); }));
}

/**
 * Specialization of `lb_constrain` to apply a container of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam L A type inheriting from `EigenBase` or a scalar.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_not_std_vector_t<L>* = nullptr>
inline auto lb_constrain(const std::vector<T>& x, const L& lb) {
  std::vector<plain_type_t<decltype(lb_constrain(x[0], lb))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], lb);
  }
  return ret;
}

/**
 * Specialization of `lb_constrain` to apply a container of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam L A type inheriting from `EigenBase` or a standard vector.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_not_std_vector_t<L>* = nullptr>
inline auto lb_constrain(const std::vector<T>& x, const L& lb,
                         return_type_t<T, L>& lp) {
  std::vector<plain_type_t<decltype(lb_constrain(x[0], lb))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], lb, lp);
  }
  return ret;
}

/**
 * Specialization of `lb_constrain` to apply a container of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam L A type inheriting from `EigenBase` or a standard vector.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L>
inline auto lb_constrain(const std::vector<T>& x, const std::vector<L>& lb) {
  check_matching_dims("lb_constrain", "x", x, "lb", lb);
  std::vector<plain_type_t<decltype(lb_constrain(x[0], lb[0]))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], lb[i]);
  }
  return ret;
}

/**
 * Specialization of `lb_constrain` to apply a container of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam L A type inheriting from `EigenBase` or a standard vector.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L>
inline auto lb_constrain(const std::vector<T>& x, const std::vector<L>& lb,
                         return_type_t<T, L>& lp) {
  check_matching_dims("lb_constrain", "x", x, "lb", lb);
  std::vector<plain_type_t<decltype(lb_constrain(x[0], lb[0]))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], lb[i], lp);
  }
  return ret;
}

/**
 * Specialization of `lb_constrain` to apply a container of lower bounds
 * elementwise to each input element. If the `Jacobian` parameter is `true`, the
 * log density accumulator is incremented with the log absolute Jacobian
 * determinant of the transform. All of the transforms are specified with their
 * Jacobians in the *Stan Reference Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam L A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in, out] lp log density accumulator
 * @return lower-bound constrained value corresponding to inputs
 */
template <bool Jacobian, typename T, typename L>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  if (Jacobian) {
    return lb_constrain(x, lb, lp);
  } else {
    return lb_constrain(x, lb);
  }
}

}  // namespace math
}  // namespace stan

#endif
