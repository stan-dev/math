#ifndef STAN_MATH_REV_CONSTRAINT_LUB_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_LUB_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/constraint/lb_constrain.hpp>
#include <stan/math/prim/constraint/lub_constrain.hpp>
#include <stan/math/rev/constraint/ub_constrain.hpp>
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
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  using std::exp;
  const auto lb_val = value_of(lb);
  const auto ub_val = value_of(ub);
  const bool is_lb_inf = lb_val == NEGATIVE_INFTY;
  const bool is_ub_inf = ub_val == INFTY;
  if (unlikely(is_ub_inf && is_lb_inf)) {
    return identity_constrain(x, ub, lb);
  } else if (unlikely(is_ub_inf)) {
    return lb_constrain(identity_constrain(x, ub), lb);
  } else if (unlikely(is_lb_inf)) {
    return ub_constrain(identity_constrain(x, lb), ub);
  } else {
    check_less("lub_constrain", "lb", lb_val, ub_val);
    auto diff = ub_val - lb_val;
    double inv_logit_x = inv_logit(value_of(x));
    return make_callback_var(
        diff * inv_logit_x + lb_val,
        [x, ub, lb, diff, inv_logit_x](auto& vi) mutable {
          if (!is_constant<T>::value) {
            forward_as<var>(x).adj()
                += vi.adj() * diff * inv_logit_x * (1.0 - inv_logit_x);
          }
          if (!is_constant<L>::value) {
            forward_as<var>(lb).adj() += vi.adj() * (1.0 - inv_logit_x);
          }
          if (!is_constant<U>::value) {
            forward_as<var>(ub).adj() += vi.adj() * inv_logit_x;
          }
        });
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
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  using std::exp;
  const auto lb_val = value_of(lb);
  const auto ub_val = value_of(ub);
  const bool is_lb_inf = lb_val == NEGATIVE_INFTY;
  const bool is_ub_inf = ub_val == INFTY;
  if (unlikely(is_ub_inf && is_lb_inf)) {
    return identity_constrain(x, ub, lb);
  } else if (unlikely(is_ub_inf)) {
    return lb_constrain(identity_constrain(x, ub), lb, lp);
  } else if (unlikely(is_lb_inf)) {
    return ub_constrain(identity_constrain(x, lb), ub, lp);
  } else {
    check_less("lub_constrain", "lb", lb_val, ub_val);
    auto neg_abs_x = -abs(value_of(x));
    auto diff = ub_val - lb_val;
    double inv_logit_x = inv_logit(value_of(x));
    lp += (log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x))));
    return make_callback_var(
        diff * inv_logit_x + lb_val,
        [x, ub, lb, diff, lp, inv_logit_x](auto& vi) mutable {
          if (!is_constant<T>::value) {
            forward_as<var>(x).adj()
                += vi.adj() * diff * inv_logit_x * (1.0 - inv_logit_x)
                   + lp.adj() * (1.0 - 2.0 * inv_logit_x);
          }
          if (!is_constant<L>::value && !is_constant<U>::value) {
            const auto one_over_diff = 1.0 / diff;
            forward_as<var>(lb).adj()
                += vi.adj() * (1.0 - inv_logit_x) + -one_over_diff * lp.adj();
            forward_as<var>(ub).adj()
                += vi.adj() * inv_logit_x + one_over_diff * lp.adj();
          } else if (!is_constant<L>::value) {
            forward_as<var>(lb).adj()
                += vi.adj() * (1.0 - inv_logit_x) + (-1.0 / diff) * lp.adj();
          } else if (!is_constant<U>::value) {
            forward_as<var>(ub).adj()
                += vi.adj() * inv_logit_x + (1.0 / diff) * lp.adj();
          }
        });
  }
}

/**
 * Specialization for Eigen matrix and scalar bounds.
 */
template <typename T, typename L, typename U, require_matrix_t<T>* = nullptr,
          require_all_stan_scalar_t<L, U>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  const auto lb_val = value_of(lb);
  const auto ub_val = value_of(ub);
  const bool is_lb_inf = lb_val == NEGATIVE_INFTY;
  const bool is_ub_inf = ub_val == INFTY;
  if (unlikely(is_ub_inf && is_lb_inf)) {
    return ret_type(identity_constrain(x, ub, lb));
  } else if (unlikely(is_ub_inf)) {
    return ret_type(lb_constrain(identity_constrain(x, ub), lb));
  } else if (unlikely(is_lb_inf)) {
    return ret_type(ub_constrain(identity_constrain(x, lb), ub));
  } else {
    arena_t<T> arena_x = x;
    check_less("lub_constrain", "lb", lb_val, ub_val);
    const auto diff = ub_val - lb_val;
    auto inv_logit_x = to_arena(inv_logit(arena_x.val().array()));
    arena_t<ret_type> ret = diff * inv_logit_x + lb_val;
    reverse_pass_callback([arena_x, ub, lb, ret, diff, inv_logit_x]() mutable {
      if (!is_constant<T>::value) {
        using T_var = arena_t<promote_scalar_t<var, T>>;
        forward_as<T_var>(arena_x).adj().array()
            += ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x);
      }
      if (!is_constant<L>::value) {
        forward_as<var>(lb).adj()
            += (ret.adj().array() * (1.0 - inv_logit_x)).sum();
      }
      if (!is_constant<U>::value) {
        forward_as<var>(ub).adj() += (ret.adj().array() * inv_logit_x).sum();
      }
    });
    return ret_type(ret);
  }
}

/**
 * Specialization for Eigen matrix and scalar bounds plus lp.
 */
template <typename T, typename L, typename U, require_matrix_t<T>* = nullptr,
          require_all_stan_scalar_t<L, U>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  const auto lb_val = value_of(lb);
  const auto ub_val = value_of(ub);
  const bool is_lb_inf = lb_val == NEGATIVE_INFTY;
  const bool is_ub_inf = ub_val == INFTY;
  if (unlikely(is_ub_inf && is_lb_inf)) {
    return ret_type(identity_constrain(x, ub, lb));
  } else if (unlikely(is_ub_inf)) {
    return ret_type(lb_constrain(identity_constrain(x, ub), lb, lp));
  } else if (unlikely(is_lb_inf)) {
    return ret_type(ub_constrain(identity_constrain(x, lb), ub, lp));
  } else {
    check_less("lub_constrain", "lb", lb_val, ub_val);
    arena_t<T> arena_x = x;
    auto neg_abs_x = to_arena(-(value_of(arena_x).array()).abs());
    auto diff = ub_val - lb_val;
    lp += (log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x)))).sum();
    auto inv_logit_x = to_arena(inv_logit(value_of(arena_x).array()));
    arena_t<ret_type> ret = diff * inv_logit_x + lb_val;
    reverse_pass_callback(
        [arena_x, ub, lb, ret, lp, diff, inv_logit_x]() mutable {
          if (!is_constant<T>::value) {
            forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x).adj().array()
                += ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x)
                   + lp.adj() * (1.0 - 2.0 * inv_logit_x);
          }
          if (!is_constant<L>::value && !is_constant<U>::value) {
            const auto lp_calc = lp.adj() * ret.size();
            const auto one_over_diff = 1.0 / diff;
            forward_as<var>(lb).adj()
                += (ret.adj().array() * (1.0 - inv_logit_x)).sum()
                   + -one_over_diff * lp_calc;
            forward_as<var>(ub).adj() += (ret.adj().array() * inv_logit_x).sum()
                                         + one_over_diff * lp_calc;
          } else if (!is_constant<L>::value) {
            forward_as<var>(lb).adj()
                += (ret.adj().array() * (1.0 - inv_logit_x)).sum()
                   + -(1.0 / diff) * lp.adj() * ret.size();
          } else if (!is_constant<U>::value) {
            forward_as<var>(ub).adj() += (ret.adj().array() * inv_logit_x).sum()
                                         + (1.0 / diff) * lp.adj() * ret.size();
          }
        });
    return ret_type(ret);
  }
}

/**
 * Specialization for Eigen matrix with matrix lower bound and scalar upper
 * bound.
 */
template <typename T, typename L, typename U,
          require_all_matrix_t<T, L>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  const auto ub_val = value_of(ub);
  const bool is_ub_inf = ub_val == INFTY;
  if (unlikely(is_ub_inf)) {
    return eval(lb_constrain(identity_constrain(x, ub), lb));
  } else {
    arena_t<T> arena_x = x;
    arena_t<L> arena_lb = lb;
    const auto lb_val = value_of(arena_lb).array().eval();
    check_less("lub_constrain", "lb", lb_val, ub_val);
    auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
    auto diff = to_arena(ub_val - lb_val);
    auto inv_logit_x = to_arena(inv_logit(value_of(arena_x).array()));
    arena_t<ret_type> ret = (is_lb_inf).select(
        ub_val - value_of(arena_x).array().exp(), diff * inv_logit_x + lb_val);
    reverse_pass_callback([arena_x, ub, arena_lb, ret, diff, inv_logit_x,
                           is_lb_inf]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using L_var = arena_t<promote_scalar_t<var, L>>;
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array() += (is_lb_inf).select(
            ret.adj().array() * -value_of(arena_x).array().exp(),
            ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x));
      }
      if (!is_constant<U>::value) {
        forward_as<var>(ub).adj()
            += (is_lb_inf)
                   .select(ret.adj().array(), ret.adj().array() * inv_logit_x)
                   .sum();
      }
      if (!is_constant<L>::value) {
        forward_as<L_var>(arena_lb).adj().array()
            += (is_lb_inf).select(0, ret.adj().array() * (1.0 - inv_logit_x));
      }
    });
    return ret_type(ret);
  }
}

/**
 * Specialization for Eigen matrix with matrix lower bound and scalar upper
 * bound plus lp.
 */
template <typename T, typename L, typename U,
          require_all_matrix_t<T, L>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          std::decay_t<return_type_t<T, L, U>>& lp) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  const auto ub_val = value_of(ub);
  const bool is_ub_inf = ub_val == INFTY;
  if (unlikely(is_ub_inf)) {
    return ret_type(lb_constrain(identity_constrain(x, ub), lb, lp));
  } else {
    arena_t<T> arena_x = x;
    arena_t<L> arena_lb = lb;
    const auto arena_x_val = to_arena(value_of(arena_x).array());
    const auto lb_val = value_of(arena_lb).array().eval();
    check_less("lub_constrain", "lb", lb_val, ub_val);
    auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
    auto diff = to_arena(ub_val - lb_val);
    auto neg_abs_x = to_arena(-arena_x_val.abs());
    auto inv_logit_x = to_arena(inv_logit(arena_x_val));
    arena_t<ret_type> ret = (is_lb_inf).select(ub_val - arena_x_val.exp(),
                                               diff * inv_logit_x + lb_val);
    lp += (is_lb_inf)
              .select(arena_x_val,
                      log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x))))
              .sum();
    reverse_pass_callback([arena_x, arena_x_val, ub, arena_lb, ret, lp, diff,
                           inv_logit_x, is_lb_inf]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using L_var = arena_t<promote_scalar_t<var, L>>;
      const auto lp_adj = lp.adj();
      if (!is_constant<T>::value) {
        const auto x_sign = arena_x_val.sign().eval();
        forward_as<T_var>(arena_x).adj().array() += (is_lb_inf).select(
            ret.adj().array() * -arena_x_val.exp() + lp_adj,
            ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x)
                + lp.adj() * (1.0 - 2.0 * inv_logit_x));
      }
      if (!is_constant<L>::value) {
        forward_as<L_var>(arena_lb).adj().array()
            += (is_lb_inf).select(0, ret.adj().array() * (1.0 - inv_logit_x)
                                         + -(1.0 / diff) * lp_adj);
      }
      if (!is_constant<U>::value) {
        forward_as<var>(ub).adj()
            += (is_lb_inf)
                   .select(ret.adj().array(), ret.adj().array() * inv_logit_x
                                                  + (1.0 / diff) * lp_adj)
                   .sum();
      }
    });
    return ret_type(ret);
  }
}

/**
 * Specialization for Eigen matrix with scalar lower bound and matrix upper
 * bound.
 */
template <typename T, typename L, typename U,
          require_all_matrix_t<T, U>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  const auto lb_val = value_of(lb);
  const bool is_lb_inf = lb_val == NEGATIVE_INFTY;
  if (unlikely(is_lb_inf)) {
    return eval(ub_constrain(identity_constrain(x, lb), ub));
  } else {
    arena_t<T> arena_x = x;
    auto arena_x_val = to_arena(value_of(arena_x));
    arena_t<U> arena_ub = ub;
    const auto ub_val = value_of(arena_ub).array().eval();
    check_less("lub_constrain", "lb", lb_val, ub_val);
    auto is_ub_inf = to_arena((ub_val == INFTY));
    auto diff = to_arena(ub_val - lb_val);
    auto inv_logit_x = to_arena(inv_logit(arena_x_val.array()));
    arena_t<ret_type> ret = (is_ub_inf).select(
        arena_x_val.array().exp() + lb_val, diff * inv_logit_x + lb_val);
    reverse_pass_callback([arena_x, arena_x_val, arena_ub, lb, ret, is_ub_inf,
                           inv_logit_x, diff]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using U_var = arena_t<promote_scalar_t<var, U>>;
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array() += (is_ub_inf).select(
            ret.adj().array() * arena_x_val.array().exp(),
            ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x));
      }
      if (!is_constant<L>::value) {
        forward_as<var>(lb).adj()
            += (is_ub_inf)
                   .select(ret.adj().array(),
                           ret.adj().array() * (1.0 - inv_logit_x))
                   .sum();
      }
      if (!is_constant<U>::value) {
        forward_as<U_var>(arena_ub).adj().array()
            += (is_ub_inf).select(0, ret.adj().array() * inv_logit_x);
      }
    });
    return ret_type(ret);
  }
}

/**
 * Specialization for Eigen matrix with scalar lower bound and matrix upper
 * bound plus lp.
 */
template <typename T, typename L, typename U,
          require_all_matrix_t<T, U>* = nullptr,
          require_stan_scalar_t<L>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          std::decay_t<return_type_t<T, L, U>>& lp) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  const auto lb_val = value_of(lb);
  const bool is_lb_inf = lb_val == NEGATIVE_INFTY;
  if (unlikely(is_lb_inf)) {
    return eval(ub_constrain(identity_constrain(x, lb), ub, lp));
  } else {
    arena_t<T> arena_x = x;
    auto arena_x_val = to_arena(value_of(arena_x));
    arena_t<U> arena_ub = ub;
    const auto& ub_val = to_ref(value_of(arena_ub));
    check_less("lub_constrain", "lb", lb_val, ub_val);
    auto is_ub_inf = to_arena((ub_val.array() == INFTY));
    auto diff = to_arena(ub_val.array() - lb_val);
    auto neg_abs_x = to_arena(-(arena_x_val.array()).abs());
    lp += (is_ub_inf)
              .select(arena_x_val.array(),
                      log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x))))
              .sum();
    auto inv_logit_x = to_arena(inv_logit(arena_x_val.array()));
    arena_t<ret_type> ret = (is_ub_inf).select(
        arena_x_val.array().exp() + lb_val, diff * inv_logit_x + lb_val);
    reverse_pass_callback([arena_x, arena_x_val, diff, inv_logit_x, arena_ub,
                           lb, ret, lp, is_ub_inf]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using U_var = arena_t<promote_scalar_t<var, U>>;
      const auto lp_adj = lp.adj();
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array() += (is_ub_inf).select(
            ret.adj().array() * arena_x_val.array().exp() + lp_adj,
            ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x)
                + lp.adj() * (1.0 - 2.0 * inv_logit_x));
      }
      if (!is_constant<L>::value) {
        forward_as<var>(lb).adj()
            += (is_ub_inf)
                   .select(ret.adj().array(),
                           ret.adj().array() * (1.0 - inv_logit_x)
                               + -(1.0 / diff) * lp_adj)
                   .sum();
      }
      if (!is_constant<U>::value) {
        forward_as<U_var>(arena_ub).adj().array() += (is_ub_inf).select(
            0, ret.adj().array() * inv_logit_x + (1.0 / diff) * lp_adj);
      }
    });
    return ret_type(ret);
  }
}

/**
 * Specialization for Eigen matrix and matrix bounds.
 */
template <typename T, typename L, typename U,
          require_all_matrix_t<T, L, U>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  arena_t<T> arena_x = x;
  auto arena_x_val = value_of(arena_x);
  arena_t<L> arena_lb = lb;
  arena_t<U> arena_ub = ub;
  auto lb_val = value_of(arena_lb).array();
  auto ub_val = value_of(arena_ub).array();
  check_less("lub_constrain", "lb", lb_val, ub_val);
  auto inv_logit_x = to_arena(inv_logit(arena_x_val.array()));
  auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
  auto is_ub_inf = to_arena((ub_val == INFTY));
  auto is_lb_ub_inf = to_arena(is_lb_inf && is_ub_inf);
  auto diff = to_arena(ub_val - lb_val);
  // if both, identity, then lb_inf -> ub_constrain, then ub_inf -> lb_constrain
  arena_t<ret_type> ret
      = (is_lb_ub_inf)
            .select(arena_x_val.array(),
                    (is_lb_inf).select(
                        ub_val - arena_x.val().array().exp(),
                        (is_ub_inf).select(arena_x_val.array().exp() + lb_val,
                                           diff * inv_logit_x + lb_val)));
  reverse_pass_callback([arena_x, arena_x_val, inv_logit_x, arena_ub, arena_lb,
                         diff, ret, is_ub_inf, is_lb_inf,
                         is_lb_ub_inf]() mutable {
    using T_var = arena_t<promote_scalar_t<var, T>>;
    using L_var = arena_t<promote_scalar_t<var, L>>;
    using U_var = arena_t<promote_scalar_t<var, U>>;
    // The most likely case is neither of them are infinity
    const bool is_none_inf = !(is_lb_inf.any() || is_ub_inf.any());
    if (is_none_inf) {
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array()
            += ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x);
      }
      if (!is_constant<L>::value) {
        forward_as<L_var>(arena_lb).adj().array()
            += ret.adj().array() * (1.0 - inv_logit_x);
      }
      if (!is_constant<U>::value) {
        forward_as<U_var>(arena_ub).adj().array()
            += ret.adj().array() * inv_logit_x;
      }
    } else {
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array()
            += (is_lb_ub_inf)
                   .select(
                       ret.adj().array(),
                       (is_lb_inf).select(
                           ret.adj().array() * -arena_x_val.array().exp(),
                           (is_ub_inf).select(
                               ret.adj().array() * arena_x_val.array().exp(),
                               ret.adj().array() * diff * inv_logit_x
                                   * (1.0 - inv_logit_x))));
      }
      if (!is_constant<U>::value) {
        forward_as<U_var>(arena_ub).adj().array() += (is_ub_inf).select(
            0, (is_lb_inf).select(ret.adj().array(),
                                  ret.adj().array() * inv_logit_x));
      }
      if (!is_constant<L>::value) {
        forward_as<L_var>(arena_lb).adj().array() += (is_lb_inf).select(
            0, (is_ub_inf).select(ret.adj().array(),
                                  ret.adj().array() * (1.0 - inv_logit_x)));
      }
    }
  });
  return ret_type(ret);
}

/**
 * Specialization for Eigen matrix and matrix bounds plus lp.
 */
template <typename T, typename L, typename U,
          require_all_matrix_t<T, L, U>* = nullptr,
          require_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub,
                          return_type_t<T, L, U>& lp) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  arena_t<T> arena_x = x;
  auto arena_x_val = value_of(arena_x);
  arena_t<L> arena_lb = lb;
  arena_t<U> arena_ub = ub;
  auto lb_val = value_of(arena_lb).array();
  auto ub_val = value_of(arena_ub).array();
  check_less("lub_constrain", "lb", lb_val, ub_val);
  auto inv_logit_x = to_arena(inv_logit(arena_x_val.array()));
  auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
  auto is_ub_inf = to_arena((ub_val == INFTY));
  auto is_lb_ub_inf = to_arena(is_lb_inf && is_ub_inf);
  auto diff = to_arena(ub_val - lb_val);
  // if both, identity, then lb_inf -> ub_constrain, then ub_inf -> lb_constrain
  arena_t<ret_type> ret
      = (is_lb_ub_inf)
            .select(arena_x_val.array(),
                    (is_lb_inf).select(
                        ub_val - arena_x.val().array().exp(),
                        (is_ub_inf).select(arena_x_val.array().exp() + lb_val,
                                           diff * inv_logit_x + lb_val)));
  auto neg_abs_x = to_arena(-(arena_x_val.array()).abs());
  lp += (is_lb_ub_inf)
            .select(
                0,
                (is_lb_inf || is_ub_inf)
                    .select(
                        arena_x_val.array(),
                        log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x)))))
            .sum();
  reverse_pass_callback([arena_x, arena_x_val, inv_logit_x, arena_ub, arena_lb,
                         diff, ret, is_ub_inf, is_lb_inf, is_lb_ub_inf,
                         lp]() mutable {
    using T_var = arena_t<promote_scalar_t<var, T>>;
    using L_var = arena_t<promote_scalar_t<var, L>>;
    using U_var = arena_t<promote_scalar_t<var, U>>;
    const auto lp_adj = lp.adj();
    // The most likely case is neither of them are infinity
    const bool is_none_inf = !(is_lb_inf.any() || is_ub_inf.any());
    if (is_none_inf) {
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array()
            += ret.adj().array() * diff * inv_logit_x * (1.0 - inv_logit_x)
               + lp.adj() * (1.0 - 2.0 * inv_logit_x);
      }
      if (!is_constant<L>::value) {
        forward_as<L_var>(arena_lb).adj().array()
            += ret.adj().array() * (1.0 - inv_logit_x) + -(1.0 / diff) * lp_adj;
      }
      if (!is_constant<U>::value) {
        forward_as<U_var>(arena_ub).adj().array()
            += ret.adj().array() * inv_logit_x + (1.0 / diff) * lp_adj;
      }
    } else {
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array()
            += (is_lb_ub_inf)
                   .select(
                       ret.adj().array(),
                       (is_lb_inf).select(
                           ret.adj().array() * -arena_x_val.array().exp()
                               + lp_adj,
                           (is_ub_inf).select(
                               ret.adj().array() * arena_x_val.array().exp()
                                   + lp_adj,
                               ret.adj().array() * diff * inv_logit_x
                                       * (1.0 - inv_logit_x)
                                   + lp.adj() * (1.0 - 2.0 * inv_logit_x))));
      }
      if (!is_constant<L>::value) {
        forward_as<L_var>(arena_lb).adj().array() += (is_lb_inf).select(
            0, (is_ub_inf).select(ret.adj().array(),
                                  ret.adj().array() * (1.0 - inv_logit_x)
                                      + -(1.0 / diff) * lp_adj));
      }
      if (!is_constant<U>::value) {
        forward_as<U_var>(arena_ub).adj().array() += (is_ub_inf).select(
            0, (is_lb_inf).select(
                   ret.adj().array(),
                   ret.adj().array() * inv_logit_x + (1.0 / diff) * lp_adj));
      }
    }
  });
  return ret_type(ret);
}

}  // namespace math
}  // namespace stan

#endif
