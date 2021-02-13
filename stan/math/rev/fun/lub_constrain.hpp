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

  // scalar

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
   template <typename T, typename L, typename U, require_all_stan_scalar_t<T, L, U>* = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
    using std::exp;
    const auto lb_val = value_of(lb);
    const auto ub_val = value_of(ub);
    const bool is_lb_inf = is_negative_infinity(lb_val);
    const bool is_ub_inf = is_positive_infinity(ub_val);
    if (unlikely(is_ub_inf && is_lb_inf)) {
      return identity_constrain(x, ub, lb);
    } else if (unlikely(is_ub_inf)) {
      return lb_constrain(identity_constrain(x, ub), lb);
    } else if (unlikely(is_lb_inf)) {
      return ub_constrain(identity_constrain(x, lb), ub);
    } else {
      check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
      return make_callback_var(add(multiply(subtract(ub_val, lb_val), inv_logit(value_of(x))), lb_val),
        [x, ub, lb](auto& vi) mutable {
          const auto exp_x = std::exp(value_of(x));
          if (!is_constant<T>::value) {
            forward_as<var>(x).adj() += vi.adj() * (exp_x / std::pow(1 + exp_x, 2)) * (value_of(ub) - value_of(lb));
          }
          if (!is_constant<L>::value) {
            forward_as<var>(lb).adj() += vi.adj() * (1 - (exp_x / (1 + exp_x)));
          }
          if (!is_constant<U>::value) {
            forward_as<var>(ub).adj() += vi.adj() * (exp_x / (1 + exp_x));
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
   template <typename T, typename L, typename U, require_all_stan_scalar_t<T, L, U>* = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub, return_type_t<T, L, U>& lp) {
    using std::exp;
    const auto lb_val = value_of(lb);
    const auto ub_val = value_of(ub);
    const bool is_lb_inf = is_negative_infinity(lb_val);
    const bool is_ub_inf = is_positive_infinity(ub_val);
    if (unlikely(is_ub_inf && is_lb_inf)) {
      return identity_constrain(x, ub, lb);
    } else if (unlikely(is_ub_inf)) {
      return lb_constrain(identity_constrain(x, ub), lb, lp);
    } else if (unlikely(is_lb_inf)) {
      return ub_constrain(identity_constrain(x, lb), ub, lp);
    } else {
      check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
      lp += add(log(subtract(ub, lb)), subtract(-abs(x), multiply(2.0, log1p_exp(-abs(x)))));
      return make_callback_var((ub_val - lb_val) * inv_logit(value_of(x)) + lb_val,
        [x, ub, lb](auto& vi) mutable {
          const auto exp_x = std::exp(value_of(x));
          if (!is_constant<T>::value) {
            forward_as<var>(x).adj() += vi.adj() * (exp_x / std::pow(1 + exp_x, 2)) * (value_of(ub) - value_of(lb));
          }
          if (!is_constant<L>::value) {
            forward_as<var>(lb).adj() += vi.adj() * (1 - (exp_x / (1 + exp_x)));
          }
          if (!is_constant<U>::value) {
            forward_as<var>(ub).adj() += vi.adj() * (exp_x / (1 + exp_x));
          }
        });
    }
  }


  /**
   * Specialization for Eigen matrix and scalar bounds.
   */
  template <typename T, typename L, typename U, require_matrix_t<T>* = nullptr,
    require_all_stan_scalar_t<L, U>* = nullptr, require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
    using std::exp;
    using ret_type = return_var_matrix_t<T, T, L, U>;
    const auto lb_val = value_of(lb);
    const auto ub_val = value_of(ub);
    const bool is_lb_inf = is_negative_infinity(lb_val);
    const bool is_ub_inf = is_positive_infinity(ub_val);
    if (unlikely(is_ub_inf && is_lb_inf)) {
      return eval(identity_constrain(x, ub, lb));
    } else if (unlikely(is_ub_inf)) {
      return eval(lb_constrain(identity_constrain(x, ub), lb));
    } else if (unlikely(is_lb_inf)) {
      return eval(ub_constrain(identity_constrain(x, lb), ub));
    } else {
      arena_t<T> arena_x = x;
      check_less("lub_constrain mat scale", "lb", lb_val, ub_val);
      arena_t<ret_type> ret = add(multiply(subtract(ub_val, lb_val), inv_logit(value_of(arena_x))), lb_val);
      reverse_pass_callback(
        [arena_x, ub, lb, ret]() mutable {
          const auto exp_x = value_of(arena_x).array().exp().eval();
          if (!is_constant<T>::value) {
            forward_as<arena_t<promote_scalar_t<var, T>>>(arena_x).adj().array() +=
            ret.adj().array() * (exp_x / (1 + exp_x).square()) * (value_of(ub) - value_of(lb));
          }
          if (!is_constant<L>::value) {
            forward_as<var>(lb).adj() += (ret.adj().array() * (1 - (exp_x / (1 + exp_x)))).sum();
          }
          if (!is_constant<U>::value) {
            forward_as<var>(ub).adj() += (ret.adj().array() * (exp_x / (1 + exp_x))).sum();
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
    const bool is_lb_inf = is_negative_infinity(lb_val);
    const bool is_ub_inf = is_positive_infinity(ub_val);
    if (unlikely(is_lb_inf)) {
      return eval(ub_constrain(x, ub));
    } else if (unlikely(is_ub_inf)) {
      return eval(lb_constrain(x, lb));
    } else {
      arena_t<T> arena_x = x;
      lp += sum(add(log(subtract(ub, lb)), subtract(-abs(x), multiply(2.0, log1p_exp(-abs(x))))));
      arena_t<ret_type> ret = (ub_val - lb_val) * inv_logit(value_of(arena_x)).array() + lb_val;
      reverse_pass_callback([ub, lb, arena_x, ret]() mutable {
          const auto exp_x = value_of(arena_x).array().exp().eval();
          auto one_plus_exp_x = (1 + exp_x).eval();
          if (!is_constant<T>::value) {
            forward_as<promote_scalar_t<var, T>>(arena_x).adj().array() += ret.adj().array() * -(exp_x / (1 + exp_x).square()) * (value_of(ub) - value_of(lb));
          }
          if (!is_constant<L>::value || !is_constant<U>::value) {
            const auto exp_over_one_plus_exp_x = (exp_x / one_plus_exp_x).eval();
            if (!is_constant<L>::value) {
              forward_as<var>(lb).adj() += (ret.adj().array() * (1 - (exp_x / (1 + exp_x)))).sum();
            }
            if (!is_constant<U>::value) {
              forward_as<var>(ub).adj() += (ret.adj().array() * (exp_x / (1 + exp_x))).sum();
            }
          }
        });
        return ret_type(ret);
    }
  }

  /**
   * Specialization for Eigen matrix with matrix lower bound and scalar upper bound.
   */
  template <typename T, typename L, typename U, require_all_matrix_t<T, L>* = nullptr,
    require_stan_scalar_t<U> * = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
    using std::exp;
    const auto ub_val = value_of(ub);
    const bool is_ub_inf = is_positive_infinity(ub_val);
    if (unlikely(is_ub_inf)) {
      return lb_constrain(x, lb);
    } else {
      arena_t<T> arena_x = x;
      arena_t<L> arena_lb = lb;
      arena_t<T> arena_ret = lub_constrain(value_of(arena_x), value_of(arena_lb), value_of(ub));
      reverse_pass_callback([ub, lb, arena_x](auto& vi) mutable {
          const auto exp_x = value_of(arena_x).array().exp().eval();
          const auto one_plus_exp_x = (1 + exp_x).eval();
          if (!is_constant<T>::value) {
            forward_as<promote_scalar_t<var, T>>(arena_x).adj().array() += vi.adj().array() * -(exp_x / (1 + exp_x).square()) * (value_of(ub) - value_of(lb));
          }
          if (!is_constant<L>::value || !is_constant<U>::value) {
            const auto one_over_one_plus_exp_x = (1 / one_plus_exp_x).eval();
            if (!is_constant<L>::value) {
              forward_as<promote_scalar_t<var, L>>(lb).adj().array() += vi.adj().array() * -one_over_one_plus_exp_x;
            }
            if (!is_constant<U>::value) {
              forward_as<var>(ub).adj() += (vi.adj().array() * one_over_one_plus_exp_x).sum();
            }
          }
        });
    }
  }


  /**
   * Specialization for Eigen matrix with matrix lower bound and scalar upper bound plus lp.
   */
  template <typename T, typename L, typename U, require_all_matrix_t<T, L>* = nullptr,
    require_stan_scalar_t<U> * = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub, std::decay_t<return_type_t<T, L, U>>& lp) {
    auto x_ref = to_ref(x);
    auto lb_ref = to_ref(lb);
    plain_type_t<T> x_ret(x.rows(), x.cols());
    for (Eigen::Index j = 0; j < x_ref.cols(); ++j) {
      for (Eigen::Index i = 0; i < x_ref.rows(); ++i) {
        x_ret.coeffRef(i, j) = lub_constrain(x_ref.coeff(i, j), lb_ref.coeff(i, j), ub, lp);
      }
    }
    return x_ret;
  }


  /**
   * Specialization for Eigen matrix with scalar lower bound and matrix upper bound.
   */
  template <typename T, typename L, typename U, require_all_matrix_t<T, U>* = nullptr,
    require_stan_scalar_t<L>* = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
    auto x_ref = to_ref(x);
    auto ub_ref = to_ref(ub);
    plain_type_t<T> x_ret(x.rows(), x.cols());
    for (Eigen::Index j = 0; j < x_ref.cols(); ++j) {
      for (Eigen::Index i = 0; i < x_ref.rows(); ++i) {
        x_ret.coeffRef(i, j) = lub_constrain(x_ref.coeff(i, j), lb, ub_ref.coeff(i, j));
      }
    }
    return x_ret;
  }


  /**
   * Specialization for Eigen matrix with scalar lower bound and matrix upper bound plus lp.
   */
  template <typename T, typename L, typename U, require_all_matrix_t<T, U>* = nullptr,
    require_stan_scalar_t<L>* = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub, std::decay_t<return_type_t<T, L, U>>& lp) {
    auto x_ref = to_ref(x);
    auto ub_ref = to_ref(ub);
    plain_type_t<T> x_ret(x.rows(), x.cols());
    for (Eigen::Index j = 0; j < x_ref.cols(); ++j) {
      for (Eigen::Index i = 0; i < x_ref.rows(); ++i) {
        x_ret.coeffRef(i, j) = lub_constrain(x_ref.coeff(i, j), lb, ub_ref.coeff(i, j), lp);
      }
    }
    return x_ret;
  }


  /**
   * Specialization for Eigen matrix and matrix bounds.
   */
  template <typename T, typename L, typename U, require_all_matrix_t<T, L, U>* = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
    auto x_ref = to_ref(x);
    auto lb_ref = to_ref(lb);
    auto ub_ref = to_ref(ub);
    plain_type_t<T> x_ret(x.rows(), x.cols());
    for (Eigen::Index j = 0; j < x_ref.cols(); ++j) {
      for (Eigen::Index i = 0; i < x_ref.rows(); ++i) {
        x_ret.coeffRef(i, j) = lub_constrain(x_ref.coeff(i, j), lb_ref.coeff(i, j), ub_ref.coeff(i, j));
      }
    }
    return x_ret;
  }


  /**
   * Specialization for Eigen matrix and matrix bounds plus lp.
   */
  template <typename T, typename L, typename U, require_all_matrix_t<T, L, U>* = nullptr,
    require_var_t<return_type_t<T, L, U>>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, const U& ub, std::decay_t<return_type_t<T, L, U>>& lp) {
    auto x_ref = to_ref(x);
    auto lb_ref = to_ref(lb);
    auto ub_ref = to_ref(ub);
    plain_type_t<T> x_ret(x.rows(), x.cols());
    for (Eigen::Index j = 0; j < x_ref.cols(); ++j) {
      for (Eigen::Index i = 0; i < x_ref.rows(); ++i) {
        x_ret.coeffRef(i, j) = lub_constrain(x_ref.coeff(i, j), lb_ref.coeff(i, j), ub_ref.coeff(i, j), lp);
      }
    }
    return x_ret;
  }


}  // namespace math
}  // namespace stan

#endif
