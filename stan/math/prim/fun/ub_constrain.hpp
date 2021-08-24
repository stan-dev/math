#ifndef STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @return matrix constrained to have upper bound
 */
template <typename T, typename U, require_all_stan_scalar_t<T, U>* = nullptr,
          require_all_not_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub) {
  if (value_of_rec(ub) == INFTY) {
    return identity_constrain(x, ub);
  } else {
    return subtract(ub, exp(x));
  }
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename T, typename U, require_all_stan_scalar_t<T, U>* = nullptr,
          require_all_not_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub,
                         std::decay_t<return_type_t<T, U>>& lp) {
  if (value_of_rec(ub) == INFTY) {
    return identity_constrain(x, ub);
  } else {
    lp += x;
    return subtract(ub, exp(x));
  }
}

/**
 * Specialization of `ub_constrain` to apply a scalar upper bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam U Scalar.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_eigen_t<T>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_all_not_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub) {
  return eval(x.unaryExpr([ub](auto&& xx) { return ub_constrain(xx, ub); }));
}

/**
 * Specialization of `ub_constrain` to apply a scalar upper bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam U Scalar.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_eigen_t<T>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_all_not_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub,
                         std::decay_t<return_type_t<T, U>>& lp) {
  return eval(
      x.unaryExpr([ub, &lp](auto&& xx) { return ub_constrain(xx, ub, lp); }));
}

/**
 * Specialization of `ub_constrain` to apply a matrix of upper bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam U A type inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_all_eigen_t<T, U>* = nullptr,
          require_all_not_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub) {
  check_matching_dims("ub_constrain", "x", x, "ub", ub);
  return eval(x.binaryExpr(
      ub, [](auto&& xx, auto&& ubb) { return ub_constrain(xx, ubb); }));
}

/**
 * Specialization of `ub_constrain` to apply a matrix of upper bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase`.
 * @tparam U A type inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_all_eigen_t<T, U>* = nullptr,
          require_all_not_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub,
                         std::decay_t<return_type_t<T, U>>& lp) {
  check_matching_dims("ub_constrain", "x", x, "ub", ub);
  return eval(x.binaryExpr(
      ub, [&lp](auto&& xx, auto&& ubb) { return ub_constrain(xx, ubb, lp); }));
}

/**
 * Specialization of `ub_constrain` to apply a scalar upper bound elementwise
 *  to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam U Scalar.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_not_std_vector_t<U>* = nullptr>
inline auto ub_constrain(const std::vector<T>& x, const U& ub) {
  std::vector<plain_type_t<decltype(ub_constrain(x[0], ub))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub);
  }
  return ret;
}

/**
 * Specialization of `ub_constrain` to apply a scalar upper bound elementwise
 *  to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam U Scalar.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_not_std_vector_t<U>* = nullptr>
inline auto ub_constrain(const std::vector<T>& x, const U& ub,
                         return_type_t<T, U>& lp) {
  std::vector<plain_type_t<decltype(ub_constrain(x[0], ub))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub, lp);
  }
  return ret;
}

/**
 * Specialization of `ub_constrain` to apply a container of upper bounds
 * elementwise to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam U A type inheriting from `EigenBase` or a standard vector.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename U>
inline auto ub_constrain(const std::vector<T>& x, const std::vector<U>& ub) {
  check_matching_dims("ub_constrain", "x", x, "ub", ub);
  std::vector<plain_type_t<decltype(ub_constrain(x[0], ub[0]))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub[i]);
  }
  return ret;
}

/**
 * Specialization of `ub_constrain` to apply a container of upper bounds
 * elementwise to each input element.
 *
 * @tparam T A Any type with a Scalar `scalar_type`.
 * @tparam U A type inheriting from `EigenBase` or a standard vector.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename U>
inline auto ub_constrain(const std::vector<T>& x, const std::vector<U>& ub,
                         return_type_t<T, U>& lp) {
  check_matching_dims("ub_constrain", "x", x, "ub", ub);
  std::vector<plain_type_t<decltype(ub_constrain(x[0], ub[0]))>> ret(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub[i], lp);
  }
  return ret;
}

/**
 * Specialization of `ub_constrain` to apply a container of upper bounds
 * elementwise to each input element. If the `Jacobian` parameter is `true`, the
 * log density accumulator is incremented with the log absolute Jacobian
 * determinant of the transform.  All of the transforms are specified with their
 * Jacobians in the *Stan Reference Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam U A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @param[in, out] lp log density accumulator
 * @return lower-bound constrained value corresponding to inputs
 */
template <bool Jacobian, typename T, typename U>
inline auto ub_constrain(const T& x, const U& ub, return_type_t<T, U>& lp) {
  if (Jacobian) {
    return ub_constrain(x, ub, lp);
  } else {
    return ub_constrain(x, ub);
  }
}

}  // namespace math
}  // namespace stan

#endif
