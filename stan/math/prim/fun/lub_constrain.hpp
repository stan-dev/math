#ifndef STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/lb_constrain.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/ub_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower and upper-bounded scalar derived by
 * transforming the specified free scalar given the specified
 * lower and upper bounds.
 *
 * <p>The transform is the transformed and scaled inverse logit,
 *
 * <p>\f$f(x) = L + (U - L) \mbox{logit}^{-1}(x)\f$
 *
 * @tparam T Scalar.
 * @tparam L Scalar.
 * @tparam U Scalar.
 * @param[in] x Free scalar to transform.
 * @param[in] lb Lower bound.
 * @param[in] ub Upper bound.
 * @return Lower- and upper-bounded scalar derived from transforming
 *   the free scalar.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U,
          require_all_stan_scalar_t<T, L, U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(T&& x, L&& lb, U&& ub) {
  const bool is_lb_inf = value_of(lb) == NEGATIVE_INFTY;
  const bool is_ub_inf = value_of(ub) == INFTY;
  if (unlikely(is_ub_inf && is_lb_inf)) {
    return identity_constrain(x, lb, ub);
  } else if (unlikely(is_ub_inf)) {
    return lb_constrain(identity_constrain(x, ub), lb);
  } else if (unlikely(is_lb_inf)) {
    return ub_constrain(identity_constrain(x, lb), ub);
  } else {
    check_less("lub_constrain", "lb", value_of(lb), value_of(ub));
    return (ub - lb) * inv_logit(x) + lb;
  }
}

/**
 * Return the lower- and upper-bounded scalar derived by
 * transforming the specified free scalar given the specified
 * lower and upper bounds and increment the specified log
 * density with the log absolute Jacobian determinant.
 *
 * <p>The transform is as defined in
 * <code>lub_constrain(T, double, double)</code>.  The log absolute
 * Jacobian determinant is given by
 *
 * <p>\f$\log \left| \frac{d}{dx} \left(
 *                L + (U-L) \mbox{logit}^{-1}(x) \right)
 *            \right|\f$
 *
 * <p>\f$ {} = \log |
 *         (U-L)
 *         \, (\mbox{logit}^{-1}(x))
 *         \, (1 - \mbox{logit}^{-1}(x)) |\f$
 *
 * <p>\f$ {} = \log (U - L) + \log (\mbox{logit}^{-1}(x))
 *                          + \log (1 - \mbox{logit}^{-1}(x))\f$
 *
 * @tparam T Scalar.
 * @tparam L Scalar.
 * @tparam U Scalar.
 * @param[in] x Free scalar to transform.
 * @param[in] lb Lower bound.
 * @param[in] ub Upper bound.
 * @param[in,out] lp Log probability scalar reference.
 * @return Lower- and upper-bounded scalar derived from transforming
 *   the free scalar.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U,
          require_all_stan_scalar_t<T, L, U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(T&& x, L&& lb, U&& ub, return_type_t<T, L, U>& lp) {
  const bool is_lb_inf = value_of(lb) == NEGATIVE_INFTY;
  const bool is_ub_inf = value_of(ub) == INFTY;
  if (unlikely(is_ub_inf && is_lb_inf)) {
    return identity_constrain(x, ub, lb);
  } else if (unlikely(is_ub_inf)) {
    return lb_constrain(identity_constrain(x, ub), lb, lp);
  } else if (unlikely(is_lb_inf)) {
    return ub_constrain(identity_constrain(x, lb), ub, lp);
  } else {
    check_less("lub_constrain", "lb", value_of(lb), value_of(ub));
    const auto diff = ub - lb;
    lp += add(log(diff), subtract(-abs(x), multiply(2.0, log1p_exp(-abs(x)))));
    return diff * inv_logit(x) + lb;
  }
}

/**
 * Overload for Eigen matrix and scalar bounds.
 */
template <typename T, typename L, typename U, require_eigen_t<T>* = nullptr,
          require_all_stan_scalar_t<L, U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  return eval(
      x.unaryExpr([ub, lb](auto&& xx) { return lub_constrain(xx, lb, ub); }));
}

/**
 * Overload for Eigen matrix and scalar bounds plus lp.
 */
template <typename T, typename L, typename U, require_eigen_t<T>* = nullptr,
          require_all_stan_scalar_t<L, U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  return eval(x.unaryExpr(
      [lb, ub, &lp](auto&& xx) { return lub_constrain(xx, lb, ub, lp); }));
}

/**
 * Overload for Eigen matrix with matrix lower bound and scalar upper
 * bound.
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, L>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  return eval(x.binaryExpr(
      lb, [ub](auto&& x, auto&& lb) { return lub_constrain(x, lb, ub); }));
}

/**
 * Overload for Eigen matrix with matrix lower bound and scalar upper
 * bound plus lp.
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, L>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  return eval(x.binaryExpr(lb, [ub, &lp](auto&& x, auto&& lb) {
    return lub_constrain(x, lb, ub, lp);
  }));
}

/**
 * Overload for Eigen matrix with scalar lower bound and matrix upper
 * bound.
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, U>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  return eval(x.binaryExpr(
      ub, [lb](auto&& x, auto&& ub) { return lub_constrain(x, lb, ub); }));
}

/**
 * Overload for Eigen matrix with scalar lower bound and matrix upper
 * bound plus lp.
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, U>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  return eval(x.binaryExpr(ub, [lb, &lp](auto&& x, auto&& ub) {
    return lub_constrain(x, lb, ub, lp);
  }));
}

/**
 * Overload for Eigen matrix and matrix bounds.
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, L, U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  auto x_ref = to_ref(x);
  auto lb_ref = to_ref(lb);
  auto ub_ref = to_ref(ub);
  promote_scalar_t<return_type_t<T, L, U>, T> x_ret(x.rows(), x.cols());
  for (Eigen::Index j = 0; j < x_ref.cols(); ++j) {
    for (Eigen::Index i = 0; i < x_ref.rows(); ++i) {
      x_ret.coeffRef(i, j) = lub_constrain(
          x_ref.coeff(i, j), lb_ref.coeff(i, j), ub_ref.coeff(i, j));
    }
  }
  return x_ret;
}

/**
 * Overload for Eigen matrix and matrix bounds plus lp.
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, L, U>* = nullptr,
          require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  auto x_ref = to_ref(x);
  auto lb_ref = to_ref(lb);
  auto ub_ref = to_ref(ub);
  promote_scalar_t<return_type_t<T, L, U>, T> x_ret(x.rows(), x.cols());
  for (Eigen::Index j = 0; j < x_ref.cols(); ++j) {
    for (Eigen::Index i = 0; i < x_ref.rows(); ++i) {
      x_ret.coeffRef(i, j) = lub_constrain(
          x_ref.coeff(i, j), lb_ref.coeff(i, j), ub_ref.coeff(i, j), lp);
    }
  }
  return x_ret;
}

/**
 * Overload for array of x and non-array lb and ub
 */
template <typename T, typename L, typename U,
          require_all_not_std_vector_t<L, U>* = nullptr>
inline auto lub_constrain(const std::vector<T>& x, const L& lb, const U& ub) {
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb, ub))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb, ub);
  }
  return ret;
}

/**
 * Overload for array of x and non-array lb and ub with lp
 */
template <typename T, typename L, typename U,
          require_all_not_std_vector_t<L, U>* = nullptr>
inline auto lub_constrain(const std::vector<T>& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb, ub))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb, ub, lp);
  }
  return ret;
}

/**
 * Overload for array of x and ub and non-array lb
 */
template <typename T, typename L, typename U,
          require_not_std_vector_t<L>* = nullptr>
inline auto lub_constrain(const std::vector<T>& x, const L& lb,
                          const std::vector<U>& ub) {
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb, ub[0]))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb, ub[i]);
  }
  return ret;
}

/**
 * Overload for array of x and ub and non-array lb with lp
 */
template <typename T, typename L, typename U,
          require_not_std_vector_t<L>* = nullptr>
inline auto lub_constrain(const std::vector<T>& x, const L& lb,
                          const std::vector<U>& ub,
                          return_type_t<T, L, U>& lp) {
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb, ub[0]))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb, ub[i], lp);
  }
  return ret;
}

/**
 * Overload for array of x and lb and non-array ub
 */
template <typename T, typename L, typename U,
          require_not_std_vector_t<U>* = nullptr>
inline auto lub_constrain(const std::vector<T>& x, const std::vector<L>& lb,
                          const U& ub) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb[0], ub))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb[i], ub);
  }
  return ret;
}

/**
 * Overload for array of x and lb and non-array ub with lp
 */
template <typename T, typename L, typename U,
          require_not_std_vector_t<U>* = nullptr>
inline auto lub_constrain(const std::vector<T>& x, const std::vector<L>& lb,
                          const U& ub, return_type_t<T, L, U>& lp) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb[0], ub))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb[i], ub, lp);
  }
  return ret;
}

/**
 * Overload for array of x, lb, and ub with lp
 */
template <typename T, typename L, typename U>
inline auto lub_constrain(const std::vector<T>& x, const std::vector<L>& lb,
                          const std::vector<U>& ub) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb[0], ub[0]))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb[i], ub[i]);
  }
  return ret;
}

/**
 * Overload for array of x, lb, and ub
 */
template <typename T, typename L, typename U>
inline auto lub_constrain(const std::vector<T>& x, const std::vector<L>& lb,
                          const std::vector<U>& ub,
                          return_type_t<T, L, U>& lp) {
  check_matching_dims("lub_constrain", "x", x, "lb", lb);
  check_matching_dims("lub_constrain", "x", x, "ub", ub);
  std::vector<plain_type_t<decltype(lub_constrain(x[0], lb[0], ub[0]))>> ret(
      x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lub_constrain(x[i], lb[i], ub[i], lp);
  }
  return ret;
}

/**
 * Return the lower and upper-bounded scalar derived by
 * transforming the specified free scalar given the specified
 * lower and upper bounds. If the `Jacobian` parameter is `true`, the log
 * density accumulator is incremented with the log absolute Jacobian determinant
 * of the transform.  All of the transforms are specified with their Jacobians
 * in the *Stan Reference Manual* chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam L A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam U A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @param[in] x Free scalar to transform
 * @param[in] lb Lower bound
 * @param[in] ub Upper bound
 * @param[in, out] lp log density accumulator
 * @return Lower- and upper-bounded scalar derived from transforming the free
 * scalar
 * @throw std::domain_error if `ub <= lb`
 */
template <bool Jacobian, typename T, typename L, typename U>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  if (Jacobian) {
    return lub_constrain(x, lb, ub, lp);
  } else {
    return lub_constrain(x, lb, ub);
  }
}

/**
 * Wrapper for tuple of bounds, simply delegates to the appropriate overload
 */
template <typename T, typename L, typename U>
inline auto lub_constrain(const T& x, const std::tuple<L, U>& bounds) {
  return lub_constrain(x, std::get<0>(bounds), std::get<1>(bounds));
}

/**
 * Wrapper for tuple of bounds, simply delegates to the appropriate overload
 */
template <typename T, typename L, typename U>
inline auto lub_constrain(const T& x, const std::tuple<L, U>& bounds,
                          return_type_t<T, L, U>& lp) {
  return lub_constrain(x, std::get<0>(bounds), std::get<1>(bounds), lp);
}

/**
 * Wrapper for tuple of bounds, simply delegates to the appropriate overload
 */
template <bool Jacobian, typename T, typename L, typename U>
inline auto lub_constrain(const T& x, const std::tuple<L, U>& bounds,
                          return_type_t<T, L, U>& lp) {
  return lub_constrain<Jacobian>(x, std::get<0>(bounds), std::get<1>(bounds),
                                 lp);
}

}  // namespace math
}  // namespace stan

#endif
