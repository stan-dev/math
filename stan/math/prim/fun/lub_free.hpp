#ifndef STAN_MATH_PRIM_FUN_LUB_FREE_HPP
#define STAN_MATH_PRIM_FUN_LUB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/divide.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <stan/math/prim/fun/lb_free.hpp>
#include <stan/math/prim/fun/ub_free.hpp>
#include <stan/math/prim/fun/subtract.hpp>

namespace stan {
namespace math {

/** lub_free
 * Return the unconstrained variable that transforms to the
 * y given the specified bounds.
 *
 * <p>The transform in `lub_constrain`,
 * is reversed by a transformed and scaled logit,
 *
 * <p>\f$f^{-1}(y) = \mbox{logit}(\frac{y - L}{U - L})\f$
 *
 * where \f$U\f$ and \f$L\f$ are the lower and upper bounds.
 *
 * @tparam T type of bounded object
 * @tparam L type of lower bound
 * @tparam U type of upper bound
 * @param y constrained value
 * @param lb lower bound
 * @param ub upper bound
 * @return the free object that transforms to the input
 *   given the bounds
 * @throw std::invalid_argument if the lower bound is greater than
 *   the upper bound, y is less than the lower bound, or y is
 *   greater than the upper bound
 */
///@{
/**
 * Overload for all scalar arguments
 */
template <typename T, typename L, typename U,
          require_not_std_vector_t<T>* = nullptr,
          require_all_stan_scalar_t<L, U>* = nullptr>
inline auto lub_free(T&& y, L&& lb, U&& ub) {
  const bool is_lb_inf = value_of(lb) == NEGATIVE_INFTY;
  const bool is_ub_inf = value_of(ub) == INFTY;
  if (unlikely(is_ub_inf && is_lb_inf)) {
    return identity_free(y, lb, ub);
  } else if (unlikely(is_ub_inf)) {
    return lb_free(identity_free(y, ub), lb);
  } else if (unlikely(is_lb_inf)) {
    return ub_free(identity_free(y, lb), ub);
  } else {
    auto&& y_ref = to_ref(std::forward<T>(y));
    check_bounded("lub_free", "Bounded variable", value_of(y_ref), value_of(lb),
                  value_of(ub));
    return eval(logit(divide(subtract(std::forward<decltype(y_ref)>(y_ref), lb),
                             subtract(ub, lb))));
  }
}

/**
 * Overload for matrix constrained variable, matrix lower bound, scalar upper
 * bound
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, L>* = nullptr,
          require_stan_scalar_t<U>* = nullptr>
inline auto lub_free(T&& y, L&& lb, U&& ub) {
  check_matching_dims("lub_free", "y", y, "lb", lb);
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& lb_ref = to_ref(std::forward<L>(lb));
  promote_scalar_t<return_type_t<T, L, U>, T> ret(y.rows(), y.cols());
  for (Eigen::Index j = 0; j < y.cols(); ++j) {
    for (Eigen::Index i = 0; i < y.rows(); ++i) {
      ret.coeffRef(i, j) = lub_free(y_ref.coeff(i, j), lb_ref.coeff(i, j), ub);
    }
  }
  return ret;
}

/**
 * Overload for matrix constrained variable, matrix upper bound, scalar lower
 * bound
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, U>* = nullptr,
          require_stan_scalar_t<L>* = nullptr>
inline auto lub_free(T&& y, L&& lb, U&& ub) {
  check_matching_dims("lub_free", "y", y, "ub", ub);
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& ub_ref = to_ref(std::forward<U>(ub));
  promote_scalar_t<return_type_t<T, L, U>, T> ret(y.rows(), y.cols());
  for (Eigen::Index j = 0; j < y.cols(); ++j) {
    for (Eigen::Index i = 0; i < y.rows(); ++i) {
      ret.coeffRef(i, j) = lub_free(y_ref.coeff(i, j), lb, ub_ref.coeff(i, j));
    }
  }
  return ret;
}

/**
 * Overload for matrix constrained variable, matrix upper bound, matrix lower
 * bound
 */
template <typename T, typename L, typename U,
          require_all_eigen_t<T, L, U>* = nullptr>
inline auto lub_free(T&& y, L&& lb, U&& ub) {
  check_matching_dims("lub_free", "y", y, "lb", lb);
  check_matching_dims("lub_free", "y", y, "ub", ub);
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& lb_ref = to_ref(std::forward<L>(lb));
  auto&& ub_ref = to_ref(std::forward<U>(ub));
  promote_scalar_t<return_type_t<T, L, U>, T> ret(y.rows(), y.cols());
  for (Eigen::Index j = 0; j < y.cols(); ++j) {
    for (Eigen::Index i = 0; i < y.rows(); ++i) {
      ret.coeffRef(i, j)
          = lub_free(y_ref.coeff(i, j), lb_ref.coeff(i, j), ub_ref.coeff(i, j));
    }
  }
  return ret;
}

/**
 * Overload for `std::vector` constrained variable
 */
template <typename T, typename L, typename U,
          require_all_not_std_vector_t<L, U>* = nullptr>
inline auto lub_free(const std::vector<T> y, const L& lb, const U& ub) {
  std::vector<decltype(lub_free(y[0], lb, ub))> ret(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    ret[i] = lub_free(y[i], lb, ub);
  }
  return ret;
}

/**
 * Overload for `std::vector` constrained variable and `std::vector` upper bound
 */
template <typename T, typename L, typename U,
          require_not_std_vector_t<L>* = nullptr>
inline auto lub_free(const std::vector<T> y, const L& lb,
                     const std::vector<U>& ub) {
  check_matching_dims("lub_free", "y", y, "ub", ub);
  std::vector<decltype(lub_free(y[0], lb, ub[0]))> ret(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    ret[i] = lub_free(y[i], lb, ub[i]);
  }
  return ret;
}

/**
 * Overload for `std::vector` constrained variable and `std::vector` lower bound
 */
template <typename T, typename L, typename U,
          require_not_std_vector_t<U>* = nullptr>
inline auto lub_free(const std::vector<T> y, const std::vector<L>& lb,
                     const U& ub) {
  check_matching_dims("lub_free", "y", y, "lb", lb);
  std::vector<decltype(lub_free(y[0], lb[0], ub))> ret(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    ret[i] = lub_free(y[i], lb[i], ub);
  }
  return ret;
}

/**
 * Overload for `std::vector` constrained variable and `std::vector` constraints
 */
template <typename T, typename L, typename U>
inline auto lub_free(const std::vector<T> y, const std::vector<L>& lb,
                     const std::vector<U>& ub) {
  check_matching_dims("lub_free", "y", y, "lb", lb);
  check_matching_dims("lub_free", "y", y, "ub", ub);
  std::vector<decltype(lub_free(y[0], lb[0], ub[0]))> ret(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    ret[i] = lub_free(y[i], lb[i], ub[i]);
  }
  return ret;
}

/**
 * Wrapper for tuple of bounds, simply delegates to the appropriate overload
 */
template <typename T, typename L, typename U>
inline auto lub_free(T&& y, const std::tuple<L, U>& bounds) {
  return lub_free(std::forward<T>(y), std::get<0>(bounds), std::get<1>(bounds));
}
///@}

}  // namespace math
}  // namespace stan
#endif
