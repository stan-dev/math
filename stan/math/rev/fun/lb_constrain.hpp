#ifndef STAN_MATH_REV_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_LB_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/fun/identity_constrain.hpp>
#include <stan/math/rev/fun/identity_free.hpp>
#include <stan/math/prim/fun/lb_constrain.hpp>
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
inline auto lb_constrain(const T& x, const L& lb) {
  using std::exp;
  const auto lb_val = value_of(lb);
  const bool is_lb_inf = is_negative_infinity(lb_val);
  if(!is_constant<T>::value && !is_constant<L>::value) {
    if (unlikely(is_lb_inf)) {
      return identity_constrain(x, lb);
    } else {
      //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
      return make_callback_var(std::exp(value_of(x)) + lb_val,
        [arena_x = var(x), arena_lb = var(lb)](auto& vi) mutable {
            arena_x.adj() += vi.adj() * std::exp(arena_x.val());
            arena_lb.adj() += vi.adj();
        });
    }
  } else if(!is_constant<T>::value) {
    if (unlikely(is_lb_inf)) {
      return identity_constrain(x, lb);
    } else {
      //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
      return make_callback_var(std::exp(value_of(x)) + lb_val,
        [arena_x = var(x)](auto& vi) mutable {
            arena_x.adj() += vi.adj() * std::exp(arena_x.val());
        });
    }
  } else {
    if (unlikely(is_lb_inf)) {
      return identity_constrain(x, lb);
    } else {
      //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
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
 * @tparam T type of Matrix
 * @tparam L type of lower bound
 * @param[in] x Unconstrained Matrix input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr,
 require_any_var_t<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, var& lp) {
  using std::exp;
  const auto lb_val = value_of(lb);
  const bool is_lb_inf = is_negative_infinity(lb_val);
  if(!is_constant<T>::value && !is_constant<L>::value) {
    if (unlikely(is_lb_inf)) {
      return identity_constrain(x, lb);
    } else {
      //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
      lp += value_of(x);
      return make_callback_var(std::exp(value_of(x)) + lb_val,
        [lp, arena_x = var(x), arena_lb = var(lb)](auto& vi) mutable {
          arena_x.adj() += vi.adj() * std::exp(arena_x.val()) + lp.adj();
          arena_lb.adj() += vi.adj();
          lp.adj() += vi.adj();
      });
    }
  } else if(!is_constant<T>::value) {
    if (unlikely(is_lb_inf)) {
      return identity_constrain(x, lb);
    } else {
      //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
      lp += value_of(x);
      return make_callback_var(std::exp(value_of(x)) + lb_val,
       [lp, arena_x = var(x)](auto& vi) mutable {
          arena_x.adj() += vi.adj() * std::exp(arena_x.val()) + lp.adj();
          lp.adj() += vi.adj();
      });
    }
  } else {
    if (unlikely(is_lb_inf)) {
      return identity_constrain(x, lb);
    } else {
      //check_less("lb_constrain", "lb", value_of(x), value_of(lb));
      lp += value_of(x);
      return make_callback_var(std::exp(value_of(x)) + lb_val,
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
             require_matrix_t<T>* = nullptr,
             require_stan_scalar_t<L>* = nullptr,
             require_any_st_var<T, L>* = nullptr>
   inline auto lb_constrain(const T& x, const L& lb) {
     using std::exp;
     using ret_type = return_var_matrix_t<T, T, L>;
     const auto lb_val = value_of(lb);
     const bool is_lb_inf = lb_val == NEGATIVE_INFTY;
     if(!is_constant<T>::value && !is_constant<L>::value) {
       if (unlikely(is_lb_inf)) {
         return identity_constrain(x, lb);
       } else {
         arena_t<promote_scalar_t<var, T>> arena_x = x;
         //check_less("lb_constrain", "lb", value_of(arena_x), value_of(lb));
         arena_t<ret_type> ret = arena_x.val().array().exp() + lb_val;
         reverse_pass_callback([arena_x, ret, arena_lb = var(lb)]() mutable {
           arena_x.adj().array() += ret.adj().array() * arena_x.val_op().array().exp();
           arena_lb.adj() += ret.adj().sum();
         });
         return ret_type(ret);
       }
     } else if(!is_constant<T>::value) {
       if (unlikely(is_lb_inf)) {
         return identity_constrain(x, lb);
       } else {
         arena_t<promote_scalar_t<var, T>> arena_x = x;
         //check_less("lb_constrain", "lb", value_of(arena_x), value_of(lb));
         arena_t<ret_type> ret = arena_x.val().array().exp() + lb_val;
         reverse_pass_callback([arena_x, ret]() mutable {
             arena_x.adj().array() += ret.adj().array() * arena_x.val_op().array().exp();
         });
         return ret_type(ret);
       }
     } else {
       if (unlikely(is_lb_inf)) {
         return identity_constrain(x, lb);
       } else {
         auto x_ref = to_ref(x);
         //check_less("lb_constrain", "lb", value_of(x_ref), value_of(lb));
         arena_t<ret_type> ret = value_of(x_ref).array().exp() + lb_val;
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
            require_matrix_t<T>* = nullptr,
            require_stan_scalar_t<L>* = nullptr,
            require_any_st_var<T, L>* = nullptr>
  inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
    using std::exp;
    using ret_type = return_var_matrix_t<T, T, L>;
    const auto lb_val = value_of(lb);
    const bool is_lb_inf = lb_val == NEGATIVE_INFTY;

    if(!is_constant<T>::value && !is_constant<L>::value) {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        arena_t<promote_scalar_t<var, T>> arena_x = x;
        //check_less("lb_constrain", "lb", value_of(arena_x), value_of(lb));
        arena_t<ret_type> ret = arena_x.val().array().exp() + lb_val;
        lp += arena_x.val().sum();
        reverse_pass_callback([arena_x, ret, lp, arena_lb = var(lb)]() mutable {
            arena_x.adj().array() += ret.adj().array() * arena_x.val().array().exp() + lp.adj();
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
        //check_less("lb_constrain", "lb", value_of(arena_x), value_of(lb));
        arena_t<ret_type> ret = arena_x.val().array().exp() + lb_val;
        lp += arena_x.val().sum();
        reverse_pass_callback([arena_x, ret, lp]() mutable {
            arena_x.adj().array() += ret.adj().array() * arena_x.val().array().exp() + lp.adj();
            lp.adj() += ret.adj().sum();
        });
        return ret_type(ret);
      }
    } else {
      if (unlikely(is_lb_inf)) {
        return identity_constrain(x, lb);
      } else {
        auto x_ref = to_ref(x);
        //check_less("lb_constrain", "lb", value_of(x_ref), value_of(lb));
        arena_t<ret_type> ret = value_of(x_ref).array().exp() + lb_val;
        lp += value_of(x).sum();
        reverse_pass_callback([ret, lp, arena_lb = var(lb)]() mutable {
            const double ret_adj = ret.adj().sum();
            arena_lb.adj() += ret_adj;
            lp.adj() += ret_adj;
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
inline auto lb_constrain(const T& x, const L& lb) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L>;
  if(!is_constant<T>::value && !is_constant<L>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    arena_t<ret_type> ret = lb_constrain(arena_x.val(), arena_lb.val());
    reverse_pass_callback([arena_x, arena_lb, ret]() mutable {
      const auto check_fin = arena_lb.val().array().isFinite().template cast<double>().eval();
      arena_x.adj().array() += ret.adj().array() * (arena_x.val().array() * check_fin).exp();
      arena_lb.adj().array() += ret.adj().array() * check_fin;
      /*
      for(size_t j = 0; j < arena_x.cols(); ++j) {
        for(size_t i = 0; i < arena_x.rows(); ++i) {
          double x_val = arena_x.val().coeff(i, j);
          double lb_val = arena_lb.val().coeff(i, j);
          double ret_adj = ret.adj().coeff(i, j);

          if (unlikely(lb_val == NEGATIVE_INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret_adj;
          } else {
            arena_x.adj().coeffRef(i, j) += ret_adj * exp(x_val);
            arena_lb.adj().coeffRef(i, j) += ret_adj;
          }
        }
      }
      */
      });

    return ret_type(ret);
  } else if(!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<double, L>> arena_lb = value_of(lb);
    arena_t<ret_type> ret = lb_constrain(arena_x.val(), arena_lb.val());

    reverse_pass_callback([arena_x, arena_lb, ret]() mutable {
      arena_x.adj().array() += ret.adj().array() * (arena_x.val().array() * arena_lb.val().array().isFinite().template cast<double>()).exp();
      /*
      for(size_t j = 0; j < arena_x.cols(); ++j) {
        for(size_t i = 0; i < arena_x.rows(); ++i) {
          double x_val = arena_x.val().coeff(i, j);
          double lb_val = arena_lb.val().coeff(i, j);
          double ret_adj = ret.adj().coeff(i, j);
          if (unlikely(lb_val == NEGATIVE_INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret_adj;
          } else {
            arena_x.adj().coeffRef(i, j) += ret_adj * exp(x_val);
          }
        }
      }
      */
    });

    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;
    arena_t<ret_type> ret = lb_constrain(value_of(x), arena_lb.val());

    reverse_pass_callback([arena_lb, ret]() mutable {
      arena_lb.adj().array() += ret.adj().array() * arena_lb.val().array().isFinite().template cast<double>();
       /*
      for(size_t j = 0; j < arena_lb.cols(); ++j) {
        for(size_t i = 0; i < arena_lb.rows(); ++i) {
          double lb_val = arena_lb.val().coeff(i, j);
          double ret_adj = ret.adj().coeff(i, j);
          if (lb_val != NEGATIVE_INFTY) {
            arena_lb.adj().coeffRef(i, j) += ret_adj;
          }
        }
      }
      */
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
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  using std::exp;
  using ret_type = return_var_matrix_t<T, T, L>;

  if(!is_constant<T>::value && !is_constant<L>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;

    double lp_val = 0.0;
    arena_t<ret_type> ret = lb_constrain(arena_x.val(), arena_lb.val(), lp_val);
    lp += lp_val;

    reverse_pass_callback([arena_x, arena_lb, ret, lp]() mutable {
      const auto check_fin = arena_lb.val().array().isFinite().template cast<double>().eval();
      arena_x.adj().array() += ret.adj().array() * (arena_x.val().array() * check_fin).exp() + lp.adj() * check_fin;
      arena_lb.adj().array() += ret.adj().array() * check_fin;
      /*
      for(size_t j = 0; j < arena_x.cols(); ++j) {
        for(size_t i = 0; i < arena_x.rows(); ++i) {
          double x_val = arena_x.val().coeff(i, j);
          double lb_val = arena_lb.val().coeff(i, j);
          double ret_adj = ret.adj().coeff(i, j);
          if (unlikely(lb_val == NEGATIVE_INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret_adj;
          } else {
            arena_x.adj().coeffRef(i, j) += ret_adj * exp(x_val) + lp_adj;
            arena_lb.adj().coeffRef(i, j) += ret_adj;
          }
        }
      }
      */
      });

    return ret_type(ret);
  } else if(!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<double, L>> arena_lb = value_of(lb);

    double lp_val = 0.0;
    arena_t<ret_type> ret = lb_constrain(arena_x.val(), arena_lb.val(), lp_val);
    lp += lp_val;

    reverse_pass_callback([arena_x, arena_lb, ret, lp]() mutable {
      double lp_adj = lp.adj();

      for(size_t j = 0; j < arena_x.cols(); ++j) {
        for(size_t i = 0; i < arena_x.rows(); ++i) {
          double x_val = arena_x.val().coeff(i, j);
          double lb_val = arena_lb.val().coeff(i, j);
          double ret_adj = ret.adj().coeff(i, j);
          if (unlikely(lb_val == NEGATIVE_INFTY)) {
            arena_x.adj().coeffRef(i, j) += ret_adj;
          } else {
            arena_x.adj().coeffRef(i, j) += ret_adj * exp(x_val) + lp_adj;
          }
        }
      }
    });

    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<var, L>> arena_lb = lb;

    double lp_val = 0.0;
    arena_t<ret_type> ret = lb_constrain(value_of(x), arena_lb.val(), lp_val);
    lp += lp_val;

    reverse_pass_callback([arena_lb, ret]() mutable {
      for(size_t j = 0; j < arena_lb.cols(); ++j) {
        for(size_t i = 0; i < arena_lb.rows(); ++i) {
          double lb_val = arena_lb.val().coeff(i, j);
          double ret_adj = ret.adj().coeff(i, j);

          if (lb_val != NEGATIVE_INFTY) {
            arena_lb.adj().coeffRef(i, j) += ret_adj;
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
