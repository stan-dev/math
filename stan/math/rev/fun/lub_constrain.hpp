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
  template <typename T1, typename T2, typename T3, require_arithmetic_t<T3>* = nullptr>
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
  template <typename T1, typename T2, typename T3, require_eigen_t<T3>* = nullptr>
  void setup_inv_logit_and_exp_x(T1& exp_x, T2& inv_logit_x, const T3& x_val) {
    for (Eigen::Index j = 0; j < x_val.cols(); ++j) {
      for (Eigen::Index i = 0; i < x_val.rows(); ++i) {
        if (x_val.coeff(i, j) < 0) {
          exp_x.coeffRef(i, j) = std::exp(x_val.coeff(i, j));
          if (x_val.coeff(i, j) > LOG_EPSILON) {
            inv_logit_x.coeffRef(i, j) = exp_x.coeff(i, j) / (1.0 + exp_x.coeff(i, j));
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
    return make_callback_var(diff * inv_logit_x + lb_val,
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
              += vi.adj() * (1.0 - inv_logit_x)
                 + -one_over_diff * lp.adj();
          forward_as<var>(ub).adj() += vi.adj() * inv_logit_x
                                       + one_over_diff * lp.adj();
        } else if (!is_constant<L>::value) {
          forward_as<var>(lb).adj()
              += vi.adj() * (1.0 - inv_logit_x)
                 + (-1.0 / diff) * lp.adj();
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
    reverse_pass_callback([arena_x, ub, lb, ret, diff, exp_x, inv_logit_x]() mutable {
      //          std::cout << "\nret adj\n" << ret.adj() << "\n";
      if (!is_constant<T>::value) {
        forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x).adj().array()
            += ret.adj().array() * (exp_x * diff) / (1.0 + exp_x).square();
        //            std::cout << "\nx adj\n" <<
        //            forward_as<arena_t<promote_scalar_t<var,
        //            T>>>(arena_x).adj() << "\n";
      }
      if (!is_constant<L>::value) {
        forward_as<var>(lb).adj()
            += (ret.adj().array() * (1.0 - inv_logit_x)).sum();
        //            std::cout << "\nlb adj\n" << forward_as<var>(lb).adj() <<
        //            "\n";
      }
      if (!is_constant<U>::value) {
        forward_as<var>(ub).adj()
            += (ret.adj().array() * inv_logit_x).sum();
        //            std::cout << "\nub adj\n" << forward_as<var>(ub).adj() <<
        //            "\n";
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
inline auto lub_constrain(const T& x, const L& lb, const U& ub, return_type_t<T, L, U>& lp) {
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
    reverse_pass_callback([arena_x, ub, lb, ret, lp, diff,
                           neg_abs_x, exp_x, inv_logit_x]() mutable {
      if (!is_constant<T>::value) {
        const auto x_sign = value_of(arena_x).array().sign().eval();
        const auto exp_neg_abs_x = (neg_abs_x).exp().eval();
        forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x).adj().array()
            += ret.adj().array() * (exp_x * diff)
                   / (1 + exp_x).square()
               + ((2.0 * exp_neg_abs_x * x_sign) / (1 + exp_neg_abs_x) - x_sign)
                     * lp.adj();
      }
      if (!is_constant<L>::value && !is_constant<U>::value) {
        const auto lp_calc = lp.adj() * ret.size();
        const auto one_over_diff = 1.0 / diff;
        forward_as<var>(lb).adj()
            += (ret.adj().array() * (1.0 - inv_logit_x)).sum() + -one_over_diff * lp_calc;
        forward_as<var>(ub).adj()
            += (ret.adj().array() * inv_logit_x).sum() + one_over_diff * lp_calc;
      } else if (!is_constant<L>::value) {
        forward_as<var>(lb).adj()
            += (ret.adj().array() * (1.0 - inv_logit_x)).sum() + -(1.0 / diff) * lp.adj() * ret.size();
      } else if (!is_constant<U>::value) {
        forward_as<var>(ub).adj()
            += (ret.adj().array() * inv_logit_x).sum() + (1.0 / diff) * lp.adj() * ret.size();
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
    arena_t<ret_type> ret = (is_lb_inf).select(
        ub_val - value_of(arena_x).array().exp(),
      diff * inv_logit(value_of(arena_x)).array() + lb_val);
    reverse_pass_callback(
        [arena_x, ub, arena_lb, ret, diff, is_lb_inf]() mutable {
          //      std::cout << "\nret adj\n" << ret.adj() << "\n";
          for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
            for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
              // In this case we do ub_constrains
              if (unlikely(is_lb_inf.coeff(i, j))) {
                const auto neg_exp_x = exp(-value_of(arena_x).coeff(i, j));
                if (!is_constant<T>::value) {
                  forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                      .adj()
                      .coeffRef(i, j)
                      += ret.adj().coeff(i, j)
                         * -std::exp(value_of(arena_x).coeff(i, j));
                }
                if (!is_constant<U>::value) {
                  forward_as<var>(ub).adj() += ret.adj().coeffRef(i, j);
                }
              } else {
                const auto neg_exp_x = std::exp(-value_of(arena_x).coeff(i, j));
                const auto one_plus_neg_exp_x = 1.0 + neg_exp_x;
                const auto ret_adj = ret.adj().coeff(i, j);
                if (!is_constant<T>::value) {
                  forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                      .adj()
                      .coeffRef(i, j)
                      += ret_adj * (neg_exp_x / std::pow(one_plus_neg_exp_x, 2))
                         * diff.coeff(i, j);
                }
                if (!is_constant<L>::value) {
                  forward_as<arena_t<promote_scalar_t<var, L>>>(arena_lb)
                      .adj()
                      .coeffRef(i, j)
                      += ret_adj * (1.0 - (1.0 / one_plus_neg_exp_x));
                }
                if (!is_constant<U>::value) {
                  forward_as<var>(ub).adj() += ret_adj * (1.0 / one_plus_neg_exp_x);
                }
              }
            }
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
inline auto lub_constrain(const T& x, const L& lb, const U& ub, std::decay_t<return_type_t<T, L, U>>& lp) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L, U>;
  const auto ub_val = value_of(ub);
  const bool is_ub_inf = ub_val == INFTY;
  if (unlikely(is_ub_inf)) {
    return eval(lb_constrain(identity_constrain(x, ub), lb, lp));
  } else {
    arena_t<T> arena_x = x;
    arena_t<L> arena_lb = lb;
    const auto lb_val = value_of(arena_lb).array().eval();
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
    auto diff = to_arena(ub_val - lb_val);
    auto neg_abs_x = to_arena(-(value_of(arena_x).array()).abs());
    arena_t<ret_type> ret = (is_lb_inf).select(
        ub_val - value_of(arena_x).array().exp(),
      diff * inv_logit(value_of(arena_x)).array() + lb_val);
   lp += (is_lb_inf).select(
     value_of(arena_x).array(),
     log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x)))).sum();
/*
    arena_t<ret_type> ret = value_of(arena_x).binaryExpr(
        value_of(arena_lb), [ub_val](const auto& x_val, const auto& lbb_val) {
          return lub_constrain(x_val, lbb_val, ub_val);
        });
        */
    reverse_pass_callback(
        [arena_x, ub, arena_lb, ret, lp, diff, neg_abs_x, is_lb_inf]() mutable {
          //      std::cout << "\nret adj\n" << ret.adj() << "\n";
          using T_var = arena_t<promote_scalar_t<var, T>>;
          using L_var = arena_t<promote_scalar_t<var, T>>;
          /*
          std::cout << "\n------------------------\n";
          std::cout << "\nx val:\n" << value_of(arena_x) <<
          "\n";
          std::cout << "\nx inv_logit:\n"
                    << inv_logit(value_of(arena_x)) << "\n";
          std::cout << "\nret adj:\n" << ret.adj() << "\n";
          std::cout << "\nlb val:\n" << value_of(arena_lb) << "\n";
          std::cout << "\nub val:\n" << value_of(ub) << "\n";
          */
          const auto lp_adj = lp.adj();
          for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
            for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
              const auto ret_adj = ret.adj().coeff(i, j);
              const auto diff_coeff = diff.coeff(i, j);
              const auto x_coeff = value_of(arena_x).coeff(i, j);
              const auto neg_exp_x_coeff = std::exp(-x_coeff);
              // In this case we do ub_constrains
              if (unlikely(is_lb_inf.coeff(i, j))) {
                if (!is_constant<T>::value) {
                  forward_as<T_var>(arena_x).adj().coeffRef(i, j) +=
                    ret_adj * -std::exp(x_coeff) + lp_adj;
                }
                if (!is_constant<U>::value) {
                  forward_as<var>(ub).adj() += ret_adj;
                  //            std::cout << "\nub adj\n" <<
                  //            forward_as<var>(ub).adj() << "\n";
                }
              } else {
                const auto one_plus_neg_exp_x_coeff = 1.0 + neg_exp_x_coeff;
                if (!is_constant<T>::value) {
                  const auto x_sign_coeff = sign(x_coeff);
                  const auto exp_neg_abs_x_coeff = std::exp(neg_abs_x.coeff(i, j));
                  forward_as<T_var>(arena_x).adj().coeffRef(i, j)
                  += ret_adj * (neg_exp_x_coeff * diff_coeff)
                         / std::pow(one_plus_neg_exp_x_coeff, 2)
                     + ((2.0 * exp_neg_abs_x_coeff * x_sign_coeff) / (1 + exp_neg_abs_x_coeff) - x_sign_coeff)
                           * lp_adj;
                }
                if (!is_constant<L>::value && !is_constant<U>::value) {
                  const auto one_over_one_plus_neg_exp_x_coeff = (1.0 / one_plus_neg_exp_x_coeff);
                  const auto one_over_diff_coeff = 1 / diff_coeff;
                  forward_as<L_var>(arena_lb).adj().coeffRef(i, j) +=
                   ret_adj * (1.0 - one_over_one_plus_neg_exp_x_coeff) +
                    -one_over_diff_coeff * lp_adj;
                  forward_as<var>(ub).adj() += ret_adj * one_over_one_plus_neg_exp_x_coeff +
                    one_over_diff_coeff * lp_adj;
                } else if (!is_constant<L>::value) {
                  forward_as<L_var>(arena_lb).adj().coeffRef(i, j) +=
                   ret_adj * (1.0 - (1.0 / one_plus_neg_exp_x_coeff)) + -(1.0 / diff_coeff) * lp_adj;
                } else if (!is_constant<U>::value) {
                  forward_as<var>(ub).adj() += ret_adj * (1.0 / one_plus_neg_exp_x_coeff) +
                  (1.0 / diff_coeff) * lp_adj;
                }
                lp.adj() += ret_adj;
              }
            }
          }
          /*
          if (!is_constant<T>::value) {
            std::cout << "\nx adj\n" <<
          forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x).adj() <<
          "\n";
          }
          if (!is_constant<L>::value) {
            std::cout << "\nlb adj\n" <<
          forward_as<arena_t<promote_scalar_t<var, L>>>(arena_lb).adj() <<
          "\n";
          }
          if (!is_constant<U>::value) {
            std::cout << "\nub adj\n" << forward_as<var>(ub).adj() <<
          "\n";
          }
          std::cout << "\nlp val:" << lp.val() << "\n";
          std::cout << "\nlp adj:" << lp.adj() << "\n";
          std::cout << "\n------------------------\n";
*/

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
    arena_t<ret_type> ret = (is_ub_inf).select(
        arena_x_val.array().exp() + lb_val,
      diff * inv_logit(arena_x_val).array() + lb_val);
    reverse_pass_callback(
        [arena_x, arena_x_val, arena_ub, lb, ret, is_ub_inf, diff]() mutable {
          using T_var = arena_t<promote_scalar_t<var, T>>;
          using U_var = arena_t<promote_scalar_t<var, U>>;
          //      std::cout << "\nret adj\n" << ret.adj() << "\n";
          for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
            for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
              const auto x_val = arena_x_val.coeff(i, j);
              const auto ret_adj = ret.adj().coeff(i, j);
              // In this case we do lb_constrains
              if (unlikely(is_ub_inf.coeff(i, j))) {
                if (!is_constant<T>::value) {
                  forward_as<T_var>(arena_x)
                      .adj()
                      .coeffRef(i, j)
                      += ret_adj
                         * std::exp(x_val);
                }
                if (!is_constant<L>::value) {
                  forward_as<var>(lb).adj() += ret.adj().coeffRef(i, j);
                  //            std::cout << "\nub adj\n" <<
                  //            forward_as<var>(ub).adj() << "\n";
                }
              } else {
                const auto neg_exp_x = std::exp(-x_val);
                const auto one_plus_neg_exp_x = 1.0 + neg_exp_x;
                const auto ret_adj = ret.adj().coeff(i, j);
                /*
                std::cout << "\nx val:\n" << value_of(arena_x).coeff(i, j) <<
                "\n"; std::cout << "\nx inv_logit:\n"
                          << inv_logit(value_of(arena_x).coeff(i, j)) << "\n";
                std::cout << "\nret adj:\n" << ret_adj << "\n";
                std::cout << "\nlb val:\n" << lb_val << "\n";
                std::cout << "\nub val:\n" << ub_val << "\n";
                */
                if (!is_constant<T>::value) {
                  forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                      .adj()
                      .coeffRef(i, j)
                      += ret_adj * (neg_exp_x / std::pow(one_plus_neg_exp_x, 2))
                         * diff.coeff(i, j);
                }
                if (!is_constant<L>::value) {
                  forward_as<var>(lb).adj() += ret_adj * (1.0 - (1.0 / one_plus_neg_exp_x));
                }
                if (!is_constant<U>::value) {
                  forward_as<arena_t<promote_scalar_t<var, U>>>(arena_ub).adj().coeffRef(i, j) += ret_adj * (1.0 / one_plus_neg_exp_x);
                }
                /*
                if (!is_constant<T>::value) {
                  std::cout << "\nx adj\n" <<
                forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x).adj() <<
                "\n";
                }
                if (!is_constant<L>::value) {
                  std::cout << "\nlb adj\n" <<
                forward_as<arena_t<promote_scalar_t<var, L>>>(arena_lb).adj() <<
                "\n";
                }
                if (!is_constant<U>::value) {
                  std::cout << "\nub adj\n" << forward_as<var>(ub).adj() <<
                "\n";
                }
                */
              }
            }
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
inline auto lub_constrain(const T& x, const L& lb, const U& ub, return_type_t<T, L, U>& lp) {
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
    const auto ub_val = value_of(arena_ub).array().eval();
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    auto is_ub_inf = to_arena((ub_val == INFTY));
    auto diff = to_arena(ub_val - lb_val);
    arena_t<ret_type> ret = (is_ub_inf).select(
        arena_x_val.array().exp() + lb_val,
      diff * inv_logit(arena_x_val).array() + lb_val);
    auto neg_abs_x = to_arena(-(arena_x_val.array()).abs());
    lp += (is_ub_inf).select(
      arena_x_val.array(),
      log(diff) + (neg_abs_x - (2.0 * log1p_exp(neg_abs_x)))).sum();
    reverse_pass_callback([arena_x, neg_abs_x, arena_x_val, diff, arena_ub, lb, ret, lp,
                           is_ub_inf]() mutable {
     using T_var = arena_t<promote_scalar_t<var, T>>;
     using U_var = arena_t<promote_scalar_t<var, U>>;
      const auto lp_adj = lp.adj();
      for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
          const auto ret_adj = ret.adj().coeff(i, j);
          const auto diff_coeff = diff.coeff(i, j);
          const auto x_coeff = arena_x_val.coeff(i, j);
          const auto neg_exp_x_coeff = std::exp(-x_coeff);
          // In this case we do lb_constrains
          if (unlikely(is_ub_inf.coeff(i, j))) {
            if (!is_constant<T>::value) {
              forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                  .adj()
                  .coeffRef(i, j)
                  += ret.adj().coeff(i, j)
                         * std::exp(value_of(arena_x).coeff(i, j))
                     + lp.adj();
              lp.adj() += ret.adj().coeff(i, j);
            }
            if (!is_constant<L>::value) {
              forward_as<var>(lb).adj() += ret.adj().coeffRef(i, j);
            }
          } else {
            const auto one_plus_neg_exp_x_coeff = 1.0 + neg_exp_x_coeff;
            if (!is_constant<T>::value) {
              const auto x_sign_coeff = sign(x_coeff);
              const auto exp_neg_abs_x_coeff = std::exp(neg_abs_x.coeff(i, j));
              forward_as<T_var>(arena_x).adj().coeffRef(i, j)
              += ret_adj * (neg_exp_x_coeff * diff_coeff)
                     / std::pow(one_plus_neg_exp_x_coeff, 2)
                 + ((2.0 * exp_neg_abs_x_coeff * x_sign_coeff) / (1 + exp_neg_abs_x_coeff) - x_sign_coeff)
                       * lp_adj;
            }
            if (!is_constant<L>::value && !is_constant<U>::value) {
              const auto one_over_one_plus_neg_exp_x_coeff = (1.0 / one_plus_neg_exp_x_coeff);
              const auto one_over_diff_coeff = 1 / diff_coeff;
              forward_as<var>(lb).adj() +=
               ret_adj * (1.0 - one_over_one_plus_neg_exp_x_coeff) +
                -one_over_diff_coeff * lp_adj;
              forward_as<U_var>(arena_ub).adj().coeffRef(i, j) += ret_adj * one_over_one_plus_neg_exp_x_coeff +
                one_over_diff_coeff * lp_adj;
            } else if (!is_constant<L>::value) {
              forward_as<var>(lb).adj() +=
               ret_adj * (1.0 - (1.0 / one_plus_neg_exp_x_coeff)) + -(1.0 / diff_coeff) * lp_adj;
            } else if (!is_constant<U>::value) {
              forward_as<U_var>(arena_ub).adj().coeffRef(i, j) += ret_adj * (1.0 / one_plus_neg_exp_x_coeff) +
              (1.0 / diff_coeff) * lp_adj;
            }
            lp.adj() += ret_adj;
          }
        }
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
  auto lb_val = value_of(arena_lb).array().eval();
  auto ub_val = value_of(arena_ub).array().eval();
  check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
  auto is_lb_inf = to_arena((lb_val == NEGATIVE_INFTY));
  auto is_ub_inf = to_arena((ub_val == INFTY));
  auto diff = to_arena(ub_val - lb_val);
  // if both, identity, then lb_inf -> ub_constrain, then ub_inf -> lb_constrain
  arena_t<ret_type> ret = (is_lb_inf && is_ub_inf).select(
    identity_constrain(arena_x_val.array(), lb, ub),
    (is_lb_inf).select(
      ub_val - arena_x.val().array().exp(),
      (is_ub_inf).select(
        arena_x_val.array().exp() + lb_val,
        diff * inv_logit(arena_x_val).array() + lb_val
      )
    )
  );
  reverse_pass_callback([arena_x, arena_x_val, arena_ub, arena_lb, diff, ret, is_ub_inf,
                         is_lb_inf]() mutable {
    using T_var = arena_t<promote_scalar_t<var, T>>;
    using L_var = arena_t<promote_scalar_t<var, L>>;
    using U_var = arena_t<promote_scalar_t<var, U>>;
    for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
      for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
        const auto x_val = arena_x_val.coeff(i, j);
        const auto ret_adj = ret.adj().coeff(i, j);
        // Need to handle both inf, inf lb, inf ub
        if (unlikely(is_ub_inf.coeff(i, j)
                     && is_lb_inf.coeff(i, j))) {
          if (!is_constant<T>::value) {
            forward_as<T_var>(arena_x)
                .adj()
                .coeffRef(i, j)
                += ret.adj().coeff(i, j);
          }
        } else if (unlikely(is_lb_inf.coeff(i, j))) {
          // In this case we do ub_constrains
          if (!is_constant<T>::value) {
            forward_as<T_var>(arena_x)
                .adj()
                .coeffRef(i, j)
                += ret.adj().coeff(i, j)
                   * -std::exp(value_of(arena_x).coeff(i, j));
          }
          if (!is_constant<U>::value) {
            // lb_constrain
            forward_as<U_var>(arena_ub).adj().coeffRef(i, j)
                += ret.adj().coeffRef(i, j);
          }
        } else if (unlikely(is_ub_inf.coeff(i, j))) {
          // In this case we do lb_constrain
          const auto neg_exp_x = exp(-value_of(arena_x).coeff(i, j));
          if (!is_constant<T>::value) {
            forward_as<T_var>(arena_x)
                .adj()
                .coeffRef(i, j)
                += ret.adj().coeff(i, j)
                   * std::exp(value_of(arena_x).coeff(i, j));
          }
          if (!is_constant<L>::value) {
            forward_as<L_var>(arena_lb).adj().coeffRef(i, j)
                += ret.adj().coeffRef(i, j);
          }
        } else {
          const auto neg_exp_x = std::exp(-x_val);
          const auto one_plus_neg_exp_x = 1.0 + neg_exp_x;
          if (!is_constant<T>::value) {
            forward_as<T_var>(arena_x)
                .adj()
                .coeffRef(i, j)
                += ret_adj * (neg_exp_x / std::pow(one_plus_neg_exp_x, 2))
                   * diff.coeff(i, j);
          }
          if (!is_constant<L>::value) {
            forward_as<L_var>(arena_lb).adj().coeffRef(i, j) += ret_adj * (1.0 - (1.0 / one_plus_neg_exp_x));
          }
          if (!is_constant<U>::value) {
            forward_as<U_var>(arena_ub).adj().coeffRef(i, j) += ret_adj * (1.0 / one_plus_neg_exp_x);
          }
        }
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
  const auto lb_val = value_of(lb);
  const auto ub_val = value_of(ub);
  const auto is_lb_inf = (lb_val.array() == NEGATIVE_INFTY).eval();
  const bool is_lb_all_inf = is_lb_inf.all();
  const auto is_ub_inf = (ub_val.array() == INFTY).eval();
  const bool is_ub_all_inf = is_ub_inf.all();
  if (unlikely(is_ub_all_inf && is_lb_all_inf)) {
    return eval(identity_constrain(x, ub, lb));
  } else if (unlikely(is_ub_all_inf)) {
    return eval(lb_constrain(identity_constrain(x, ub), lb, lp));
  } else if (unlikely(is_lb_all_inf)) {
    return eval(ub_constrain(identity_constrain(x, lb), ub, lp));
  } else {
    check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
    arena_t<T> arena_x = x;
    arena_t<L> arena_lb = lb;
    arena_t<U> arena_ub = ub;
    arena_t<decltype(is_lb_inf)> arena_is_lb_inf = is_lb_inf;
    arena_t<decltype(is_ub_inf)> arena_is_ub_inf = is_ub_inf;
    Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime> ret_val(
        x.rows(), x.cols());
    double lp_val = 0.0;
    for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
      for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
        ret_val.coeffRef(i, j) = lub_constrain(
            value_of(arena_x).coeff(i, j), value_of(arena_lb).coeff(i, j),
            value_of(arena_ub).coeff(i, j), lp_val);
      }
    }
    lp += lp_val;
    arena_t<ret_type> ret = ret_val;
    reverse_pass_callback([arena_x, arena_ub, arena_lb, ret, arena_is_ub_inf,
                           arena_is_lb_inf, lp]() mutable {
      for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
          // Need to handle both inf, inf lb, inf ub
          if (unlikely(arena_is_ub_inf.coeff(i, j)
                       && arena_is_lb_inf.coeff(i, j))) {
            if (!is_constant<T>::value) {
              forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                  .adj()
                  .coeffRef(i, j)
                  += ret.adj().coeff(i, j);
            }
          } else if (unlikely(arena_is_lb_inf.coeff(i, j))) {
            // In this case we do ub_constrains
            const auto neg_exp_x = exp(-value_of(arena_x).coeff(i, j));
            if (!is_constant<T>::value) {
              forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                  .adj()
                  .coeffRef(i, j)
                  += ret.adj().coeff(i, j)
                         * -std::exp(value_of(arena_x).coeff(i, j))
                     + lp.adj();
              lp.adj() += ret.adj().coeff(i, j);
            }
            if (!is_constant<U>::value) {
              forward_as<arena_t<promote_scalar_t<var, U>>>(arena_ub)
                  .adj()
                  .coeffRef(i, j)
                  += ret.adj().coeffRef(i, j);
            }
          } else if (unlikely(arena_is_ub_inf.coeff(i, j))) {
            // In this case we do lb_constrain
            const auto neg_exp_x = exp(-value_of(arena_x).coeff(i, j));
            if (!is_constant<T>::value) {
              forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                  .adj()
                  .coeffRef(i, j)
                  += ret.adj().coeff(i, j)
                         * std::exp(value_of(arena_x).coeff(i, j))
                     + lp.adj();
              lp.adj() += ret.adj().coeff(i, j);
            }
            if (!is_constant<L>::value) {
              forward_as<arena_t<promote_scalar_t<var, L>>>(arena_lb)
                  .adj()
                  .coeffRef(i, j)
                  += ret.adj().coeffRef(i, j);
            }
          } else {
            const auto neg_exp_x = exp(-value_of(arena_x).coeff(i, j));
            const auto lb_val = value_of(arena_lb).coeff(i, j);
            const auto ub_val = value_of(arena_ub).coeff(i, j);
            const auto ret_adj = ret.adj().coeff(i, j);
            const auto x_val = value_of(arena_x).coeff(i, j);
            const auto exp_x = std::exp(value_of(arena_x).coeff(i, j));
            if (!is_constant<T>::value) {
              forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x)
                  .adj()
                  .coeffRef(i, j)
                  += ret_adj * (exp_x / std::pow(1 + exp_x, 2))
                         * (ub_val - lb_val)
                     + (((2.0 * exp(-abs(x_val)) * sign(x_val))
                             / (1 + exp(-abs(x_val)))
                         - sign(x_val))
                        * lp.adj());
            }
            if (!is_constant<L>::value) {
              forward_as<arena_t<promote_scalar_t<var, U>>>(arena_lb)
                  .adj()
                  .coeffRef(i, j)
                  += ret_adj * (1 - (exp_x / (1 + exp_x)))
                     + ((-1.0 / (ub_val - lb_val)) * lp.adj());
            }
            if (!is_constant<U>::value) {
              forward_as<arena_t<promote_scalar_t<var, U>>>(arena_ub)
                  .adj()
                  .coeffRef(i, j)
                  += ret_adj * (exp_x / (1 + exp_x))
                     + (1.0 / (ub_val - lb_val) * lp.adj());
            }
          }
        }
      }
    });
    return ret_type(ret);
  }
}

}  // namespace math
}  // namespace stan

#endif
