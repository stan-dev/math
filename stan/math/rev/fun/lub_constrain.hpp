#ifndef STAN_MATH_REV_FUN_LUB_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_LUB_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/fun/lb_constrain.hpp>
#include <stan/math/prim/fun/lub_constrain.hpp>
#include <stan/math/rev/fun/ub_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
/**
 * For scalar input values calculate the underflow adjusted
 * `exp(-x)` and inverse logit.
 * @tparam T1 Scalar
 * @tparam T3 Scalar
 * @tparam T3 Scalar
 * @param[out] exp_x Returns the value of underflow corrected `exp(-x)`
 * @param[out] inv_logit_x Returns the value of underflow corrected inverse
 * logit.
 * @param[in] x_val The values of `x` to calculate both of the above from.
 */
template <typename T1, typename T2, typename T3,
          require_arithmetic_t<T3>* = nullptr>
void setup_inv_logit_and_exp_x(T1& exp_x, T2& inv_logit_x, const T3& x_val) {
  if (x_val < 0) {
    exp_x = std::exp(x_val);
    if (x_val > LOG_EPSILON) {
      inv_logit_x = exp_x / (1.0 + exp_x);
    } else {
      inv_logit_x = exp_x;
    }
  } else {
    exp_x = std::exp(-x_val);
    inv_logit_x = 1.0 / (1.0 + exp_x);
  }
}
/**
 * For Matrix input values calculate the underflow adjusted
 * `exp(-x)` and inverse logit.
 * @tparam T1 Eigen Array
 * @tparam T3 Eigen Array
 * @tparam T3 Eigen Array
 * @param[out] exp_x Returns the value of underflow corrected `exp(-x)`
 * @param[out] inv_logit_x Returns the value of underflow corrected inverse
 * logit.
 * @param[in] x_val The values of `x` to calculate both of the above from.
 */
template <typename T1, typename T2, typename T3, require_eigen_t<T3>* = nullptr>
void setup_inv_logit_and_exp_x(T1& exp_x, T2& inv_logit_x, const T3& x_val) {
  for (Eigen::Index j = 0; j < x_val.cols(); ++j) {
    for (Eigen::Index i = 0; i < x_val.rows(); ++i) {
      if (x_val.coeff(i, j) < 0) {
        exp_x.coeffRef(i, j) = std::exp(x_val.coeff(i, j));
        if (x_val.coeff(i, j) > LOG_EPSILON) {
          inv_logit_x.coeffRef(i, j)
              = exp_x.coeff(i, j) / (1.0 + exp_x.coeff(i, j));
        } else {
          inv_logit_x.coeffRef(i, j) = exp_x.coeff(i, j);
        }
      } else {
        exp_x.coeffRef(i, j) = std::exp(-x_val.coeff(i, j));
        inv_logit_x.coeffRef(i, j) = 1.0 / (1.0 + exp_x.coeff(i, j));
      }
    }
  }
}
}  // namespace internal

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
 * @tparam T type of Matrix
 * @tparam L type of lower bound
 * @param[in] x Unconstrained Matrix input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
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
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto diff = ub_val - lb_val;
    double x_val = value_of(x);
    double exp_x = 0.0;
    double inv_logit_x = 0.0;
    internal::setup_inv_logit_and_exp_x(exp_x, inv_logit_x, x_val);
    return make_callback_var(
        diff * inv_logit_x + lb_val,
        [x, ub, lb, diff, exp_x, inv_logit_x](auto& vi) mutable {
          if (!is_constant<T>::value) {
            forward_as<var>(x).adj()
                += vi.adj() * (exp_x * diff) / std::pow(1.0 + exp_x, 2);
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
 * @tparam T type of Matrix
 * @tparam L type of lower bound
 * @param[in] x Unconstrained Matrix input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
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
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto neg_abs_x = -abs(value_of(x));
    auto diff = ub_val - lb_val;
    lp += (log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x))));
    double x_val = value_of(x);
    double exp_x = 0.0;
    double inv_logit_x = 0.0;
    internal::setup_inv_logit_and_exp_x(exp_x, inv_logit_x, x_val);
    return make_callback_var(
        diff * inv_logit_x + lb_val,
        [x, ub, lb, diff, neg_abs_x, lp, exp_x, inv_logit_x](auto& vi) mutable {
          if (!is_constant<T>::value) {
            const auto x_sign = sign(x);
            const auto exp_neg_abs_x = exp(neg_abs_x);
            forward_as<var>(x).adj()
                += vi.adj() * ((exp_x * diff) / std::pow(1.0 + exp_x, 2))
                   + ((2.0 * exp_neg_abs_x * x_sign) / (1.0 + exp_neg_abs_x)
                      - x_sign)
                         * lp.adj();
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
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    const auto diff = ub_val - lb_val;
    auto x_val = value_of(x).array().eval();
    arena_t<plain_type_t<decltype(x_val)>> exp_x(x.rows(), x.cols());
    arena_t<plain_type_t<decltype(x_val)>> inv_logit_x(x.rows(), x.cols());
    internal::setup_inv_logit_and_exp_x(exp_x, inv_logit_x, x_val.array());
    arena_t<ret_type> ret = diff * inv_logit_x + lb_val;
    reverse_pass_callback([arena_x, ub, lb, ret, diff, exp_x,
                           inv_logit_x]() mutable {
      if (!is_constant<T>::value) {
        using T_var = arena_t<promote_scalar_t<var, T>>;
        forward_as<T_var>(arena_x).adj().array()
            += ret.adj().array() * (exp_x * diff) / (1.0 + exp_x).square();
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
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    arena_t<T> arena_x = x;
    auto neg_abs_x = to_arena(-(value_of(arena_x).array()).abs());
    auto diff = ub_val - lb_val;
    const auto& x_val = to_ref(value_of(x).array());
    arena_t<plain_type_t<decltype(x_val)>> exp_x(x.rows(), x.cols());
    arena_t<plain_type_t<decltype(x_val)>> inv_logit_x(x.rows(), x.cols());
    internal::setup_inv_logit_and_exp_x(exp_x, inv_logit_x, x_val.val());
    lp += (log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x)))).sum();
    arena_t<ret_type> ret = diff * inv_logit_x + lb_val;
    reverse_pass_callback([arena_x, ub, lb, ret, lp, diff, neg_abs_x, exp_x,
                           inv_logit_x]() mutable {
      if (!is_constant<T>::value) {
        const auto x_sign = value_of(arena_x).array().sign().eval();
        const auto exp_neg_abs_x = (neg_abs_x).exp().eval();
        forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x).adj().array()
            += ret.adj().array() * (exp_x * diff) / (1 + exp_x).square()
               + ((2.0 * exp_neg_abs_x * x_sign) / (1 + exp_neg_abs_x) - x_sign)
                     * lp.adj();
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
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
    auto diff = to_arena(ub_val - lb_val);
    const auto& x_val = to_ref(value_of(x).array());
    arena_t<plain_type_t<decltype(x_val)>> exp_neg_x(x.rows(), x.cols());
    arena_t<plain_type_t<decltype(x_val)>> inv_logit_x(x.rows(), x.cols());
    internal::setup_inv_logit_and_exp_x(exp_neg_x, inv_logit_x, x_val.val());
    arena_t<ret_type> ret = (is_lb_inf).select(
        ub_val - value_of(arena_x).array().exp(), diff * inv_logit_x + lb_val);
    reverse_pass_callback([arena_x, ub, arena_lb, ret, diff, inv_logit_x,
                           exp_neg_x, is_lb_inf]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using L_var = arena_t<promote_scalar_t<var, L>>;
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array() += (is_lb_inf).select(
            ret.adj().array() * -value_of(arena_x).array().exp(),
            ret.adj().array()
                * ((exp_neg_x * diff) / (1 + exp_neg_x).square()));
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
    const auto x_val = value_of(arena_x).array();
    const auto lb_val = value_of(arena_lb).array().eval();
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
    auto diff = to_arena(ub_val - lb_val);
    auto neg_abs_x = to_arena(-(value_of(arena_x).array()).abs());
    arena_t<plain_type_t<decltype(x_val)>> exp_neg_x(x.rows(), x.cols());
    arena_t<plain_type_t<decltype(x_val)>> inv_logit_x(x.rows(), x.cols());
    internal::setup_inv_logit_and_exp_x(exp_neg_x, inv_logit_x, x_val);
    arena_t<ret_type> ret = (is_lb_inf).select(
        ub_val - value_of(arena_x).array().exp(), diff * inv_logit_x + lb_val);
    lp += (is_lb_inf)
              .select(value_of(arena_x).array(),
                      log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x))))
              .sum();
    reverse_pass_callback([arena_x, ub, arena_lb, ret, lp, diff, neg_abs_x,
                           exp_neg_x, inv_logit_x, is_lb_inf]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using L_var = arena_t<promote_scalar_t<var, L>>;
      const auto lp_adj = lp.adj();
      if (!is_constant<T>::value) {
        const auto x_sign = value_of(arena_x).array().sign().eval();
        const auto exp_neg_abs_x = neg_abs_x.exp().eval();
        forward_as<T_var>(arena_x).adj().array() += (is_lb_inf).select(
            ret.adj().array() * -value_of(arena_x).array().exp() + lp_adj,
            ret.adj().array() * ((exp_neg_x * diff) / (1 + exp_neg_x).square())
                + ((2.0 * exp_neg_abs_x * x_sign) / (1 + exp_neg_abs_x)
                   - x_sign)
                      * lp_adj);
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
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto is_ub_inf = to_arena((ub_val == INFTY));
    auto diff = to_arena(ub_val - lb_val);
    arena_t<plain_type_t<decltype(arena_x_val.array())>> exp_neg_x(x.rows(),
                                                                   x.cols());
    arena_t<plain_type_t<decltype(arena_x_val.array())>> inv_logit_x(x.rows(),
                                                                     x.cols());
    internal::setup_inv_logit_and_exp_x(exp_neg_x, inv_logit_x,
                                        arena_x_val.array());
    arena_t<ret_type> ret = (is_ub_inf).select(
        arena_x_val.array().exp() + lb_val, diff * inv_logit_x + lb_val);
    reverse_pass_callback([arena_x, arena_x_val, arena_ub, lb, ret, is_ub_inf,
                           exp_neg_x, inv_logit_x, diff]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using U_var = arena_t<promote_scalar_t<var, U>>;
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array() += (is_ub_inf).select(
            ret.adj().array() * arena_x_val.array().exp(),
            ret.adj().array()
                * ((exp_neg_x * diff) / (1 + exp_neg_x).square()));
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
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto is_ub_inf = to_arena((ub_val.array() == INFTY));
    auto diff = to_arena(ub_val.array() - lb_val);
    auto neg_abs_x = to_arena(-(arena_x_val.array()).abs());
    lp += (is_ub_inf)
              .select(arena_x_val.array(),
                      log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x))))
              .sum();
    arena_t<plain_type_t<decltype(arena_x_val.array())>> exp_neg_x(x.rows(),
                                                                   x.cols());
    arena_t<plain_type_t<decltype(arena_x_val.array())>> inv_logit_x(x.rows(),
                                                                     x.cols());
    internal::setup_inv_logit_and_exp_x(exp_neg_x, inv_logit_x,
                                        arena_x_val.array());
    arena_t<ret_type> ret = (is_ub_inf).select(
        arena_x_val.array().exp() + lb_val, diff * inv_logit_x + lb_val);
    reverse_pass_callback([arena_x, neg_abs_x, arena_x_val, diff, exp_neg_x,
                           inv_logit_x, arena_ub, lb, ret, lp,
                           is_ub_inf]() mutable {
      using T_var = arena_t<promote_scalar_t<var, T>>;
      using U_var = arena_t<promote_scalar_t<var, U>>;
      const auto lp_adj = lp.adj();
      if (!is_constant<T>::value) {
        const auto x_sign = arena_x_val.array().sign().eval();
        const auto exp_neg_abs_x = neg_abs_x.exp().eval();
        forward_as<T_var>(arena_x).adj().array() += (is_ub_inf).select(
            ret.adj().array() * arena_x_val.array().exp() + lp_adj,
            ret.adj().array() * (exp_neg_x * diff) / (1 + exp_neg_x).square()
                + ((2.0 * exp_neg_abs_x * x_sign) / (1 + exp_neg_abs_x)
                   - x_sign)
                      * lp_adj);
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
  check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
  using plain_x_array = plain_type_t<decltype(arena_x_val.array())>;
  arena_t<plain_x_array> exp_neg_x(x.rows(), x.cols());
  arena_t<plain_x_array> inv_logit_x(x.rows(), x.cols());
  internal::setup_inv_logit_and_exp_x(exp_neg_x, inv_logit_x,
                                      arena_x_val.array());
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
  reverse_pass_callback([arena_x, arena_x_val, exp_neg_x, inv_logit_x, arena_ub,
                         arena_lb, diff, ret, is_ub_inf, is_lb_inf,
                         is_lb_ub_inf]() mutable {
    using T_var = arena_t<promote_scalar_t<var, T>>;
    using L_var = arena_t<promote_scalar_t<var, L>>;
    using U_var = arena_t<promote_scalar_t<var, U>>;
    // The most likely case is neither of them are infinity
    const bool is_none_inf = !(is_lb_inf.any() || is_ub_inf.any());
    if (is_none_inf) {
      if (!is_constant<T>::value) {
        forward_as<T_var>(arena_x).adj().array() += ret.adj().array()
                                                    * (exp_neg_x * diff)
                                                    / (1 + exp_neg_x).square();
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
                               ret.adj().array() * (exp_neg_x * diff)
                                   / (1 + exp_neg_x).square())));
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
  check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
  using plain_x_array = plain_type_t<decltype(arena_x_val.array())>;
  arena_t<plain_x_array> exp_neg_x(x.rows(), x.cols());
  arena_t<plain_x_array> inv_logit_x(x.rows(), x.cols());
  internal::setup_inv_logit_and_exp_x(exp_neg_x, inv_logit_x,
                                      arena_x_val.array());
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
  reverse_pass_callback([arena_x, arena_x_val, exp_neg_x, inv_logit_x,
                         neg_abs_x, arena_ub, arena_lb, diff, ret, is_ub_inf,
                         is_lb_inf, is_lb_ub_inf, lp]() mutable {
    using T_var = arena_t<promote_scalar_t<var, T>>;
    using L_var = arena_t<promote_scalar_t<var, L>>;
    using U_var = arena_t<promote_scalar_t<var, U>>;
    const auto lp_adj = lp.adj();
    // The most likely case is neither of them are infinity
    const bool is_none_inf = !(is_lb_inf.any() || is_ub_inf.any());
    if (is_none_inf) {
      if (!is_constant<T>::value) {
        const auto x_sign = arena_x_val.array().sign().eval();
        const auto exp_neg_abs_x = neg_abs_x.exp().eval();
        forward_as<T_var>(arena_x).adj().array()
            += ret.adj().array() * (exp_neg_x * diff) / (1 + exp_neg_x).square()
               + ((2.0 * exp_neg_abs_x * x_sign) / (1 + exp_neg_abs_x) - x_sign)
                     * lp_adj;
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
        const auto x_sign = arena_x_val.array().sign().eval();
        const auto exp_neg_abs_x = neg_abs_x.exp().eval();
        forward_as<T_var>(arena_x).adj().array()
            += (is_lb_ub_inf)
                   .select(ret.adj().array(),
                           (is_lb_inf).select(
                               ret.adj().array() * -arena_x_val.array().exp()
                                   + lp_adj,
                               (is_ub_inf).select(
                                   ret.adj().array() * arena_x_val.array().exp()
                                       + lp_adj,
                                   ret.adj().array() * (exp_neg_x * diff)
                                           / (1 + exp_neg_x).square()
                                       + ((2.0 * exp_neg_abs_x * x_sign)
                                              / (1 + exp_neg_abs_x)
                                          - x_sign)
                                             * lp_adj)));
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
