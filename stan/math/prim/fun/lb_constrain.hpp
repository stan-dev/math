#ifndef STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
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
 * @tparam T A scalar or matrix.
 * @tparam L A scalar or matrix.
 * @param[in] x Unconstrained input
 * @param[in] lb Lower bound
 * @return Constrained matrix
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L> * = nullptr, require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb) {

  if (unlikely(is_negative_infinity(lb))) {
    return identity_constrain(x, lb);
  } else {
    //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
    return add(exp(x), lb);
  }
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * @tparam T A scalar or matrix.
 * @tparam L A scalar or matrix.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L> * = nullptr, require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  if (is_negative_infinity(lb)) {
    return identity_constrain(x, lb);
  } else {
    //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
    lp += x;
    return add(exp(x), lb);
  }
}

template <typename T, typename L, require_eigen_t<T>* = nullptr,
  require_stan_scalar_t<L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& ub) {
  return x.unaryExpr([ub](auto&& xx) {
    return lb_constrain(xx, ub);
  });
}

template <typename T, typename L, require_all_eigen_t<T, L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& ub) {
  return x.binaryExpr(ub, [](auto&& xx, auto&& ubb) {
    return lb_constrain(xx, ubb);
  });
}

template <typename T, typename L, require_eigen_t<T>* = nullptr,
  require_stan_scalar_t<L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  return eval(x.unaryExpr([ub, &lp](auto&& xx) {
    return lb_constrain(xx, ub, lp);
  }));
}

template <typename T, typename L, require_all_eigen_t<T, L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  return eval(x.binaryExpr(ub, [&lp](auto&& xx, auto&& ubb) {
    return lb_constrain(xx, ubb, lp);
  }));
}

// VEC

template <typename T, typename L,
  require_stan_scalar_t<L>* = nullptr>
inline auto lb_constrain(const std::vector<T>& x, const L& ub) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], ub);
  }
  return ret;
}

template <typename T, typename L, require_container_t<L>* = nullptr>
inline auto lb_constrain(const std::vector<T>& x, const L& ub) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], ub[i]);
  }
  return ret;
}

template <typename T, typename L,
  require_stan_scalar_t<L>* = nullptr>
inline auto lb_constrain(const std::vector<T>& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], ub, lp);
  }
  return ret;
}

template <typename T, typename L, require_container_t<L>* = nullptr>
inline auto lb_constrain(const std::vector<T>& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = lb_constrain(x[i], ub[i], lp);
  }
  return ret;
}


}  // namespace math
}  // namespace stan

#endif
