#ifndef STAN_MATH_REV_FUN_lub_constrain_HPP
#define STAN_MATH_REV_FUN_lub_constrain_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/lub_constrain.hpp>
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
   template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr,
    require_any_var_t<T, L>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb) {
    using std::exp;
    const auto lb_val = value_of(lb);
    const bool is_lb_inf = is_positive_infinity(lb_val);
    if(!is_constant<T>::value && !is_constant<L>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        return make_callback_var(lb_val - std::exp(value_of(x)),
          [arena_x = var(x), arena_lb = var(lb)](auto& vi) mutable {
              arena_x.adj() += vi.adj() * -std::exp(arena_x.val());
              arena_lb.adj() += vi.adj();
          });
      }
    } else if(!is_constant<T>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        return make_callback_var(lb_val - std::exp(value_of(x)),
          [arena_x = var(x)](auto& vi) mutable {
              arena_x.adj() += vi.adj() * -std::exp(arena_x.val());
          });
      }
    } else {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        return make_callback_var(lb_val - std::exp(value_of(x)),
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
   * @tparam T type of Matrix
   * @tparam L type of lower bound
   * @param[in] x Unconstrained Matrix input
   * @param[in] lb lower bound on constrained output
   * @return lower bound constrained value corresponding to inputs
   */
  template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr,
   require_any_var_t<T, L>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
    using std::exp;
    const auto lb_val = value_of(lb);
    const bool is_lb_inf = lb_val == INFTY;
    if(!is_constant<T>::value && !is_constant<L>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        lp += value_of(x);
        return make_callback_var(value_of(lb) - std::exp(value_of(x)),
          [lp, arena_x = var(x), arena_lb = var(lb)](auto& vi) mutable {
            arena_x.adj() += vi.adj() * -std::exp(arena_x.val()) + lp.adj();
            arena_lb.adj() += vi.adj();
            lp.adj() += vi.adj();
        });
      }
    } else if(!is_constant<T>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        lp += value_of(x);
        return make_callback_var(value_of(lb) - std::exp(value_of(x)),
          [lp, arena_x = var(x)](auto& vi) mutable {
            arena_x.adj() += vi.adj() * -std::exp(arena_x.val()) + lp.adj();
            lp.adj() += vi.adj();
        });
      }
    } else {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        lp += value_of(x);
        return make_callback_var(value_of(lb) - std::exp(value_of(x)),
          [lp, arena_lb = var(lb)](auto& vi) mutable {
            arena_lb.adj() += vi.adj();
            lp.adj() += vi.adj();
        });
      }
    }
  }



// matrix and scalar
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
  template <typename T, typename L,
            require_any_matrix_t<T, L>* = nullptr,
            require_any_st_var<T, L>* = nullptr,
            require_stan_scalar_t<L>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb) {
    using std::exp;
    using ret_type = return_var_matrix_t<T, T, L>;
    const auto lb_val = value_of(lb);
    const bool is_lb_inf = lb_val == INFTY;
    if(!is_constant<T>::value && !is_constant<L>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        arena_t<promote_scalar_t<var, T>> arena_x = x;
        arena_t<ret_type> ret = lb_val - arena_x.val().array().exp();
        reverse_pass_callback([arena_x, ret, arena_lb = var(lb)]() mutable {
          arena_x.adj().array() += ret.adj().array() * -arena_x.val_op().array().exp();
          arena_lb.adj() += ret.adj().sum();
        });
        return ret_type(ret);
      }
    } else if(!is_constant<T>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        arena_t<promote_scalar_t<var, T>> arena_x = x;
        arena_t<ret_type> ret = lb_val - arena_x.val().array().exp();
        reverse_pass_callback([arena_x, ret]() mutable {
            arena_x.adj().array() += ret.adj().array() * -arena_x.val_op().array().exp();
        });
        return ret_type(ret);
      }
    } else {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        arena_t<ret_type> ret = lb_val - value_of(x).array().exp();
        reverse_pass_callback([ret, arena_lb = var(lb)]() mutable {
            arena_lb.adj() += ret.adj().sum();
        });
        return ret_type(ret);
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
  template <typename T, typename L,
            require_any_matrix_t<T, L>* = nullptr,
            require_any_st_var<T, L>* = nullptr,
            require_stan_scalar_t<L>* = nullptr>
  inline auto lub_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
    using std::exp;
    using ret_type = return_var_matrix_t<T, T, L>;
    const auto lb_val = value_of(lb);
    const bool is_lb_inf = lb_val == INFTY;
    if(!is_constant<T>::value && !is_constant<L>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        arena_t<promote_scalar_t<var, T>> arena_x = x;
        arena_t<promote_scalar_t<var, L>> arena_lb = lb;
        double lp_val = 0.0;
        arena_t<ret_type> ret = lb_val - arena_x.val().array().exp();
        lp += arena_x.val().sum();
        reverse_pass_callback([arena_x, arena_lb, ret, lp]() mutable {
            arena_x.adj().array() += ret.adj().array() * -arena_x.val().array().exp() + lp.adj();
            const double ret_adj_sum = ret.adj().sum();
            arena_lb.adj() += ret_adj_sum;
            lp.adj() += ret_adj_sum;
        });
        return ret_type(ret);
      }
    } else if(!is_constant<T>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        arena_t<promote_scalar_t<var, T>> arena_x = x;
        arena_t<promote_scalar_t<double, L>> arena_lb = value_of(lb);
        double lp_val = 0.0;
        arena_t<ret_type> ret = lb_val - arena_x.val().array().exp();
        lp += arena_x.val().sum();
        reverse_pass_callback([arena_x, arena_lb, ret, lp]() mutable {
            arena_x.adj().array() += ret.adj().array() * -arena_x.val().array().exp() + lp.adj();
            lp.adj() += ret.adj().sum();
        });
        return ret_type(ret);
      }
    } else {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        arena_t<promote_scalar_t<var, L>> arena_lb = lb;
        arena_t<ret_type> ret = lb_val - value_of(x).array().exp();
        lp += value_of(x).sum();
        reverse_pass_callback([arena_lb, ret, lp]() mutable {
          const double ret_adj_sum = ret.adj().sum();
          arena_lb.adj() += ret_adj_sum;
          lp.adj() += ret_adj_sum;
        });
        return ret_type(ret);
      }
    }
  }

// MATRIX

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
template <typename T, typename L,
          require_all_matrix_t<T, L>* = nullptr,
          require_any_st_var<T, L>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L>;

  if(!is_constant<T>::value && !is_constant<L>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    arena_t<ret_type> ret = lub_constrain(arena_x.val(), arena_lb.val());
    reverse_pass_callback([arena_x, arena_lb, ret]() mutable {
      for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
          if (unlikely(arena_lb.val().coeff(i, j) == INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret.adj().coeff(i, j);
          } else {
            const double ret_adj = ret.adj().coeff(i, j);
            arena_x.adj().coeffRef(i, j) += ret_adj * -std::exp(arena_x.val().coeff(i, j));
            arena_lb.adj().coeffRef(i, j) += ret_adj;
          }
        }
      }
    });
    return ret_type(ret);
  } else if(!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<double, L>> arena_lb = value_of(lb);
    arena_t<ret_type> ret = lub_constrain(arena_x.val(), arena_lb.val());
    reverse_pass_callback([arena_x, arena_lb, ret]() mutable {
      for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
          if (unlikely(arena_lb.val().coeff(i, j) == INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret.adj().coeff(i, j);
          } else {
            arena_x.adj().coeffRef(i, j) += ret.adj().coeff(i, j) * -std::exp(arena_x.val().coeff(i, j));
          }
        }
      }
    });
    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    arena_t<ret_type> ret = lub_constrain(value_of(x), arena_lb.val());
    reverse_pass_callback([arena_lb, ret]() mutable {
      for (Eigen::Index j = 0; j < arena_lb.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_lb.rows(); ++i) {
          if (arena_lb.val().coeff(i, j) != INFTY) {
            arena_lb.adj().coeffRef(i, j) += ret.adj().coeff(i, j);
          }
        }
      }
    });

    return ret_type(ret);
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
template <typename T, typename L,
          require_all_matrix_t<T, L>* = nullptr,
          require_any_st_var<T, L>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L>;
  if(!is_constant<T>::value && !is_constant<L>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    double lp_val = 0.0;
    arena_t<ret_type> ret = lub_constrain(arena_x.val(), arena_lb.val(), lp_val);
    lp += lp_val;
    reverse_pass_callback([arena_x, arena_lb, ret, lp]() mutable {
      double lp_adj = lp.adj();
      for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
          if (unlikely(arena_lb.val().coeff(i, j) == INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret.adj().coeff(i, j);
          } else {
            const double ret_adj = ret.adj().coeff(i, j);
            arena_x.adj().coeffRef(i, j) += ret_adj * -exp(arena_x.val().coeff(i, j)) + lp_adj;
            arena_lb.adj().coeffRef(i, j) += ret_adj;
          }
        }
      }
    });
    return ret_type(ret);
  } else if(!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<double, L>> arena_lb = value_of(lb);
    double lp_val = 0.0;
    arena_t<ret_type> ret = lub_constrain(arena_x.val(), arena_lb.val(), lp_val);
    lp += lp_val;
    reverse_pass_callback([arena_x, arena_lb, ret, lp]() mutable {
      double lp_adj = lp.adj();
      for (Eigen::Index j = 0; j < arena_x.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_x.rows(); ++i) {
          if (unlikely(arena_lb.val().coeff(i, j) == INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret.adj().coeff(i, j);
          } else {
            arena_x.adj().coeffRef(i, j) += ret.adj().coeff(i, j) * -std::exp(arena_x.val().coeff(i, j)) + lp_adj;
          }
        }
      }
    });
    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    double lp_val = 0.0;
    arena_t<ret_type> ret = lub_constrain(value_of(x), arena_lb.val(), lp_val);
    lp += lp_val;
    reverse_pass_callback([arena_lb, ret]() mutable {
      for (Eigen::Index j = 0; j < arena_lb.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_lb.rows(); ++i) {
          if (unlikely(arena_lb.val().coeff(i, j) != INFTY)) {
            arena_lb.adj().coeffRef(i, j) += ret.adj().coeff(i, j);
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
