#ifndef STAN_MATH_REV_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_LB_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/lb_constrain.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/rev/fun/identity_constrain.hpp>
#include <stan/math/rev/fun/identity_free.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>
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
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T Scalar
 * @tparam L Scalar
 * @param[in] x Unconstrained input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr,
          require_any_var_t<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb) {
  const auto lb_val = value_of(lb);
  if (unlikely(lb_val == NEGATIVE_INFTY)) {
    return identity_constrain(x, lb);
  } else {
    if (!is_constant<T>::value && !is_constant<L>::value) {
      auto exp_x = std::exp(value_of(x));
      return make_callback_var(
          exp_x + lb_val,
          [arena_x = var(x), arena_lb = var(lb), exp_x](auto& vi) mutable {
            arena_x.adj() += vi.adj() * exp_x;
            arena_lb.adj() += vi.adj();
          });
    } else if (!is_constant<T>::value) {
      auto exp_x = std::exp(value_of(x));
      return make_callback_var(exp_x + lb_val,
                               [arena_x = var(x), exp_x](auto& vi) mutable {
                                 arena_x.adj() += vi.adj() * exp_x;
                               });
    } else {
      return make_callback_var(std::exp(value_of(x)) + lb_val,
                               [arena_lb = var(lb)](auto& vi) mutable {
                                 arena_lb.adj() += vi.adj();
                               });
    }
  }
}

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
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T Scalar
 * @tparam L Scalar
 * @param[in] x Unconstrained input
 * @param[in] lb lower bound on constrained output
 * @param[in,out] lp reference to log probability to increment
 * @return lower bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr,
          require_any_var_t<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, var& lp) {
  const auto lb_val = value_of(lb);
  if (unlikely(lb_val == NEGATIVE_INFTY)) {
    return identity_constrain(x, lb);
  } else {
    lp += value_of(x);
    if (!is_constant<T>::value && !is_constant<L>::value) {
      auto exp_x = std::exp(value_of(x));
      return make_callback_var(
          exp_x + lb_val,
          [lp, arena_x = var(x), arena_lb = var(lb), exp_x](auto& vi) mutable {
            arena_x.adj() += vi.adj() * exp_x + lp.adj();
            arena_lb.adj() += vi.adj();
          });
    } else if (!is_constant<T>::value) {
      auto exp_x = std::exp(value_of(x));
      return make_callback_var(exp_x + lb_val,
                               [lp, arena_x = var(x), exp_x](auto& vi) mutable {
                                 arena_x.adj() += vi.adj() * exp_x + lp.adj();
                               });
    } else {
      return make_callback_var(std::exp(value_of(x)) + lb_val,
                               [lp, arena_lb = var(lb)](auto& vi) mutable {
                                 arena_lb.adj() += vi.adj();
                               });
    }
  }
}

/**
 * Specialization of `lb_constrain` to apply a scalar lower bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam L Scalar.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_matrix_t<T>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_any_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb) {
  using ret_type = return_var_matrix_t<T, T, L>;
  const auto lb_val = value_of(lb);
  if (unlikely(lb_val == NEGATIVE_INFTY)) {
    return ret_type(identity_constrain(x, lb));
  } else {
    if (!is_constant<T>::value && !is_constant<L>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto exp_x = to_arena(arena_x.val().array().exp());
      arena_t<ret_type> ret = exp_x + lb_val;
      reverse_pass_callback(
          [arena_x, ret, exp_x, arena_lb = var(lb)]() mutable {
            arena_x.adj().array() += ret.adj().array() * exp_x;
            arena_lb.adj() += ret.adj().sum();
          });
      return ret_type(ret);
    } else if (!is_constant<T>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto exp_x = to_arena(arena_x.val().array().exp());
      arena_t<ret_type> ret = exp_x + lb_val;
      reverse_pass_callback([arena_x, ret, exp_x]() mutable {
        arena_x.adj().array() += ret.adj().array() * exp_x;
      });
      return ret_type(ret);
    } else {
      arena_t<ret_type> ret = value_of(x).array().exp() + lb_val;
      reverse_pass_callback([ret, arena_lb = var(lb)]() mutable {
        arena_lb.adj() += ret.adj().sum();
      });
      return ret_type(ret);
    }
  }
}

/**
 * Specialization of `lb_constrain` to apply a scalar lower bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam L Scalar.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_matrix_t<T>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_any_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  using ret_type = return_var_matrix_t<T, T, L>;
  const auto lb_val = value_of(lb);
  if (unlikely(lb_val == NEGATIVE_INFTY)) {
    return ret_type(identity_constrain(x, lb));
  } else {
    if (!is_constant<T>::value && !is_constant<L>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto exp_x = to_arena(arena_x.val().array().exp());
      arena_t<ret_type> ret = exp_x + lb_val;
      lp += arena_x.val().sum();
      reverse_pass_callback(
          [arena_x, ret, lp, arena_lb = var(lb), exp_x]() mutable {
            arena_x.adj().array() += ret.adj().array() * exp_x + lp.adj();
            arena_lb.adj() += ret.adj().sum();
          });
      return ret_type(ret);
    } else if (!is_constant<T>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto exp_x = to_arena(arena_x.val().array().exp());
      arena_t<ret_type> ret = exp_x + lb_val;
      lp += arena_x.val().sum();
      reverse_pass_callback([arena_x, ret, exp_x, lp]() mutable {
        arena_x.adj().array() += ret.adj().array() * exp_x + lp.adj();
      });
      return ret_type(ret);
    } else {
      const auto& x_ref = to_ref(x);
      lp += sum(x_ref);
      arena_t<ret_type> ret = value_of(x_ref).array().exp() + lb_val;
      reverse_pass_callback([ret, lp, arena_lb = var(lb)]() mutable {
        arena_lb.adj() += ret.adj().sum();
      });
      return ret_type(ret);
    }
  }
}

/**
 * Specialization of `lb_constrain` to apply a matrix of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam L A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_matrix_t<T, L>* = nullptr,
          require_any_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb) {
  check_matching_dims("lb_constrain", "x", x, "lb", lb);
  using ret_type = return_var_matrix_t<T, T, L>;
  if (!is_constant<T>::value && !is_constant<L>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    auto is_not_inf_lb = to_arena((arena_lb.val().array() != NEGATIVE_INFTY));
    auto precomp_x_exp = to_arena((arena_x.val().array()).exp());
    arena_t<ret_type> ret = (is_not_inf_lb)
                                .select(precomp_x_exp + arena_lb.val().array(),
                                        arena_x.val().array());
    reverse_pass_callback([arena_x, arena_lb, ret, is_not_inf_lb,
                           precomp_x_exp]() mutable {
      arena_x.adj().array()
          += (is_not_inf_lb)
                 .select(ret.adj().array() * precomp_x_exp, ret.adj().array());
      arena_lb.adj().array() += (is_not_inf_lb).select(ret.adj().array(), 0);
    });
    return ret_type(ret);
  } else if (!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    auto lb_ref = to_ref(value_of(lb));
    auto is_not_inf_lb = to_arena((lb_ref.array() != NEGATIVE_INFTY));
    auto precomp_x_exp = to_arena((arena_x.val().array()).exp());
    arena_t<ret_type> ret
        = (is_not_inf_lb)
              .select(precomp_x_exp + lb_ref.array(), arena_x.val().array());
    reverse_pass_callback([arena_x, ret, is_not_inf_lb,
                           precomp_x_exp]() mutable {
      arena_x.adj().array()
          += (is_not_inf_lb)
                 .select(ret.adj().array() * precomp_x_exp, ret.adj().array());
    });
    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    const auto x_ref = to_ref(value_of(x));
    auto is_not_inf_lb = to_arena((arena_lb.val().array() != NEGATIVE_INFTY));
    arena_t<ret_type> ret
        = (is_not_inf_lb)
              .select(x_ref.array().exp() + arena_lb.val().array(),
                      x_ref.array());
    reverse_pass_callback([arena_lb, ret, is_not_inf_lb]() mutable {
      arena_lb.adj().array() += (is_not_inf_lb).select(ret.adj().array(), 0);
    });
    return ret_type(ret);
  }
}

/**
 * Specialization of `lb_constrain` to apply a matrix of lower bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam L A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_matrix_t<T, L>* = nullptr,
          require_any_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  check_matching_dims("lb_constrain", "x", x, "lb", lb);
  using ret_type = return_var_matrix_t<T, T, L>;
  if (!is_constant<T>::value && !is_constant<L>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    auto is_not_inf_lb = to_arena((arena_lb.val().array() != NEGATIVE_INFTY));
    auto exp_x = to_arena(arena_x.val().array().exp());
    arena_t<ret_type> ret
        = (is_not_inf_lb)
              .select(exp_x + arena_lb.val().array(), arena_x.val().array());
    lp += (is_not_inf_lb).select(arena_x.val(), 0).sum();
    reverse_pass_callback(
        [arena_x, arena_lb, ret, lp, exp_x, is_not_inf_lb]() mutable {
          const auto lp_adj = lp.adj();
          for (size_t j = 0; j < arena_x.cols(); ++j) {
            for (size_t i = 0; i < arena_x.rows(); ++i) {
              double ret_adj = ret.adj().coeff(i, j);
              if (likely(is_not_inf_lb.coeff(i, j))) {
                arena_x.adj().coeffRef(i, j)
                    += ret_adj * exp_x.coeff(i, j) + lp_adj;
                arena_lb.adj().coeffRef(i, j) += ret_adj;
              } else {
                arena_x.adj().coeffRef(i, j) += ret_adj;
              }
            }
          }
        });
    return ret_type(ret);
  } else if (!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    auto lb_val = value_of(lb).array();
    auto is_not_inf_lb = to_arena((lb_val != NEGATIVE_INFTY));
    auto exp_x = to_arena(arena_x.val().array().exp());
    arena_t<ret_type> ret
        = (is_not_inf_lb).select(exp_x + lb_val, arena_x.val().array());
    auto lp_old = lp;
    lp += (is_not_inf_lb).select(arena_x.val(), 0).sum();
    reverse_pass_callback([arena_x, ret, exp_x, lp, is_not_inf_lb]() mutable {
      const auto lp_adj = lp.adj();
      for (size_t j = 0; j < arena_x.cols(); ++j) {
        for (size_t i = 0; i < arena_x.rows(); ++i) {
          if (likely(is_not_inf_lb.coeff(i, j))) {
            const double ret_adj = ret.adj().coeff(i, j);
            arena_x.adj().coeffRef(i, j)
                += ret_adj * exp_x.coeff(i, j) + lp_adj;
          } else {
            arena_x.adj().coeffRef(i, j) += ret.adj().coeff(i, j);
          }
        }
      }
    });
    return ret_type(ret);
  } else {
    auto x_val = to_ref(value_of(x)).array();
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    auto is_not_inf_lb = to_arena((arena_lb.val().array() != NEGATIVE_INFTY));
    arena_t<ret_type> ret
        = (is_not_inf_lb).select(x_val.exp() + arena_lb.val().array(), x_val);
    lp += (is_not_inf_lb).select(x_val, 0).sum();
    reverse_pass_callback([arena_lb, ret, is_not_inf_lb]() mutable {
      arena_lb.adj().array()
          += ret.adj().array() * is_not_inf_lb.template cast<double>();
    });

    return ret_type(ret);
  }
}

}  // namespace math
}  // namespace stan

#endif
