#ifndef STAN_MATH_REV_FUN_UB_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_UB_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/ub_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * @tparam T Scalar
 * @tparam U Scalar
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @return matrix constrained to have upper bound
 */
template <typename T, typename U, require_all_stan_scalar_t<T, U>* = nullptr,
          require_any_var_t<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub) {
  const auto ub_val = value_of(ub);
  if (unlikely(ub_val == INFTY)) {
    return identity_constrain(x, ub);
  } else {
    if (!is_constant<T>::value && !is_constant<U>::value) {
      auto neg_exp_x = -std::exp(value_of(x));
      return make_callback_var(
          ub_val + neg_exp_x,
          [arena_x = var(x), arena_ub = var(ub), neg_exp_x](auto& vi) mutable {
            const auto vi_adj = vi.adj();
            arena_x.adj() += vi_adj * neg_exp_x;
            arena_ub.adj() += vi_adj;
          });
    } else if (!is_constant<T>::value) {
      auto neg_exp_x = -std::exp(value_of(x));
      return make_callback_var(ub_val + neg_exp_x,
                               [arena_x = var(x), neg_exp_x](auto& vi) mutable {
                                 arena_x.adj() += vi.adj() * neg_exp_x;
                               });
    } else {
      return make_callback_var(ub_val - std::exp(value_of(x)),
                               [arena_ub = var(ub)](auto& vi) mutable {
                                 arena_ub.adj() += vi.adj();
                               });
    }
  }
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in, out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename T, typename U, require_all_stan_scalar_t<T, U>* = nullptr,
          require_any_var_t<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub, return_type_t<T, U>& lp) {
  const auto ub_val = value_of(ub);
  const bool is_ub_inf = ub_val == INFTY;
  if (!is_constant<T>::value && !is_constant<U>::value) {
    if (unlikely(is_ub_inf)) {
      return identity_constrain(x, ub);
    } else {
      lp += value_of(x);
      auto neg_exp_x = -std::exp(value_of(x));
      return make_callback_var(value_of(ub) + neg_exp_x,
                               [lp, arena_x = var(x), arena_ub = var(ub),
                                neg_exp_x](auto& vi) mutable {
                                 const auto vi_adj = vi.adj();
                                 arena_x.adj() += vi_adj * neg_exp_x + lp.adj();
                                 arena_ub.adj() += vi_adj;
                               });
    }
  } else if (!is_constant<T>::value) {
    if (unlikely(is_ub_inf)) {
      return identity_constrain(x, ub);
    } else {
      lp += value_of(x);
      auto neg_exp_x = -std::exp(value_of(x));
      return make_callback_var(
          value_of(ub) + neg_exp_x,
          [lp, arena_x = var(x), neg_exp_x](auto& vi) mutable {
            arena_x.adj() += vi.adj() * neg_exp_x + lp.adj();
          });
    }
  } else {
    if (unlikely(is_ub_inf)) {
      return identity_constrain(x, ub);
    } else {
      lp += value_of(x);
      return make_callback_var(value_of(ub) - std::exp(value_of(x)),
                               [arena_ub = var(ub)](auto& vi) mutable {
                                 arena_ub.adj() += vi.adj();
                               });
    }
  }
}

/**
 * Specialization of `ub_constrain` to apply a scalar upper bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam U Scalar.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_matrix_t<T>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_any_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub) {
  using ret_type = return_var_matrix_t<T, T, U>;
  const auto ub_val = value_of(ub);
  if (unlikely(ub_val == INFTY)) {
    return ret_type(identity_constrain(x, ub));
  } else {
    if (!is_constant<T>::value && !is_constant<U>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto arena_neg_exp_x = to_arena(-arena_x.val().array().exp());
      arena_t<ret_type> ret = ub_val + arena_neg_exp_x;
      reverse_pass_callback(
          [arena_x, arena_neg_exp_x, ret, arena_ub = var(ub)]() mutable {
            arena_x.adj().array() += ret.adj().array() * arena_neg_exp_x;
            arena_ub.adj() += ret.adj().sum();
          });
      return ret_type(ret);
    } else if (!is_constant<T>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto arena_neg_exp_x = to_arena(-arena_x.val().array().exp());
      arena_t<ret_type> ret = ub_val + arena_neg_exp_x;
      reverse_pass_callback([arena_x, arena_neg_exp_x, ret]() mutable {
        arena_x.adj().array() += ret.adj().array() * arena_neg_exp_x;
      });
      return ret_type(ret);
    } else {
      arena_t<ret_type> ret = ub_val - value_of(x).array().exp();
      reverse_pass_callback([ret, arena_ub = var(ub)]() mutable {
        arena_ub.adj() += ret.adj().sum();
      });
      return ret_type(ret);
    }
  }
}

/**
 * Specialization of `ub_constrain` to apply a scalar upper bound elementwise
 *  to each input.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam U Scalar.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_matrix_t<T>* = nullptr,
          require_stan_scalar_t<U>* = nullptr,
          require_any_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub, return_type_t<T, U>& lp) {
  using ret_type = return_var_matrix_t<T, T, U>;
  const auto ub_val = value_of(ub);
  if (unlikely(ub_val == INFTY)) {
    return ret_type(identity_constrain(x, ub));
  } else {
    if (!is_constant<T>::value && !is_constant<U>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto arena_neg_exp_x = to_arena(-arena_x.val().array().exp());
      arena_t<ret_type> ret = ub_val + arena_neg_exp_x;
      lp += arena_x.val().sum();
      reverse_pass_callback([arena_x, arena_neg_exp_x, ret, lp,
                             arena_ub = var(ub)]() mutable {
        arena_x.adj().array() += ret.adj().array() * arena_neg_exp_x + lp.adj();
        arena_ub.adj() += ret.adj().sum();
      });
      return ret_type(ret);
    } else if (!is_constant<T>::value) {
      arena_t<promote_scalar_t<var, T>> arena_x = x;
      auto arena_neg_exp_x = to_arena(-arena_x.val().array().exp());
      arena_t<ret_type> ret = ub_val + arena_neg_exp_x;
      lp += arena_x.val().sum();
      reverse_pass_callback([arena_x, arena_neg_exp_x, ret, lp]() mutable {
        arena_x.adj().array() += ret.adj().array() * arena_neg_exp_x + lp.adj();
      });
      return ret_type(ret);
    } else {
      auto x_ref = to_ref(value_of(x));
      arena_t<ret_type> ret = ub_val - x_ref.array().exp();
      lp += x_ref.sum();
      reverse_pass_callback([ret, arena_ub = var(ub)]() mutable {
        arena_ub.adj() += ret.adj().sum();
      });
      return ret_type(ret);
    }
  }
}

/**
 * Specialization of `ub_constrain` to apply a matrix of upper bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam U A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_all_matrix_t<T, U>* = nullptr,
          require_any_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub) {
  check_matching_dims("ub_constrain", "x", x, "ub", ub);
  using ret_type = return_var_matrix_t<T, T, U>;
  if (!is_constant<T>::value && !is_constant<U>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, U>> arena_ub = ub;
    auto ub_val = to_ref(arena_ub.val());
    auto is_not_inf_ub = to_arena((ub_val.array() != INFTY));
    auto neg_exp_x = to_arena(-arena_x.val().array().exp());
    arena_t<ret_type> ret
        = (is_not_inf_ub)
              .select(ub_val.array() + neg_exp_x, arena_x.val().array());
    reverse_pass_callback([arena_x, neg_exp_x, arena_ub, ret,
                           is_not_inf_ub]() mutable {
      arena_x.adj().array()
          += (is_not_inf_ub)
                 .select(ret.adj().array() * neg_exp_x, ret.adj().array());
      arena_ub.adj().array() += (is_not_inf_ub).select(ret.adj().array(), 0.0);
    });
    return ret_type(ret);
  } else if (!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    auto ub_val = to_ref(value_of(ub));
    auto is_not_inf_ub = to_arena((ub_val.array() != INFTY));
    auto neg_exp_x = to_arena(-arena_x.val().array().exp());
    arena_t<ret_type> ret
        = (is_not_inf_ub)
              .select(ub_val.array() + neg_exp_x, arena_x.val().array());
    reverse_pass_callback([arena_x, neg_exp_x, ret, is_not_inf_ub]() mutable {
      arena_x.adj().array()
          += (is_not_inf_ub)
                 .select(ret.adj().array() * neg_exp_x, ret.adj().array());
    });
    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<var, U>> arena_ub = to_arena(ub);
    auto is_not_inf_ub
        = to_arena((arena_ub.val().array() != INFTY).template cast<double>());
    auto&& x_ref = to_ref(value_of(x).array());
    arena_t<ret_type> ret
        = (is_not_inf_ub).select(arena_ub.val().array() - x_ref.exp(), x_ref);
    reverse_pass_callback([arena_ub, ret, is_not_inf_ub]() mutable {
      arena_ub.adj().array() += ret.adj().array() * is_not_inf_ub;
    });
    return ret_type(ret);
  }
}

/**
 * Specialization of `ub_constrain` to apply a matrix of upper bounds
 * elementwise to each input element.
 *
 * @tparam T A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @tparam U A type inheriting from `EigenBase` or a `var_value` with inner type
 * inheriting from `EigenBase`.
 * @param[in] x unconstrained input
 * @param[in] ub upper bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return upper-bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_all_matrix_t<T, U>* = nullptr,
          require_any_st_var<T, U>* = nullptr>
inline auto ub_constrain(const T& x, const U& ub, return_type_t<T, U>& lp) {
  check_matching_dims("ub_constrain", "x", x, "ub", ub);
  using ret_type = return_var_matrix_t<T, T, U>;
  if (!is_constant<T>::value && !is_constant<U>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    arena_t<promote_scalar_t<var, U>> arena_ub = ub;
    auto ub_val = to_ref(arena_ub.val());
    auto is_not_inf_ub = to_arena((ub_val.array() != INFTY));
    auto neg_exp_x = to_arena(-arena_x.val().array().exp());
    arena_t<ret_type> ret
        = (is_not_inf_ub)
              .select(ub_val.array() + neg_exp_x, arena_x.val().array());
    lp += (is_not_inf_ub).select(arena_x.val().array(), 0).sum();
    reverse_pass_callback([arena_x, neg_exp_x, arena_ub, ret, lp,
                           is_not_inf_ub]() mutable {
      arena_x.adj().array()
          += (is_not_inf_ub)
                 .select(ret.adj().array() * neg_exp_x + lp.adj(),
                         ret.adj().array());
      arena_ub.adj().array() += (is_not_inf_ub).select(ret.adj().array(), 0.0);
    });
    return ret_type(ret);
  } else if (!is_constant<T>::value) {
    arena_t<promote_scalar_t<var, T>> arena_x = x;
    auto ub_val = to_ref(value_of(ub));
    auto is_not_inf_ub = to_arena((ub_val.array() != INFTY));
    auto neg_exp_x = to_arena(-arena_x.val().array().exp());
    arena_t<ret_type> ret
        = (is_not_inf_ub)
              .select(ub_val.array() + neg_exp_x, arena_x.val().array());
    lp += (is_not_inf_ub).select(arena_x.val().array(), 0).sum();
    reverse_pass_callback(
        [arena_x, neg_exp_x, ret, lp, is_not_inf_ub]() mutable {
          arena_x.adj().array()
              += (is_not_inf_ub)
                     .select(ret.adj().array() * neg_exp_x + lp.adj(),
                             ret.adj().array());
        });
    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<var, U>> arena_ub = to_arena(ub);
    auto is_not_inf_ub
        = to_arena((arena_ub.val().array() != INFTY).template cast<double>());
    auto&& x_ref = to_ref(value_of(x).array());
    arena_t<ret_type> ret
        = (is_not_inf_ub).select(arena_ub.val().array() - x_ref.exp(), x_ref);
    lp += (is_not_inf_ub).select(x_ref, 0).sum();
    reverse_pass_callback([arena_ub, ret, is_not_inf_ub]() mutable {
      arena_ub.adj().array() += ret.adj().array() * is_not_inf_ub;
    });
    return ret_type(ret);
  }
}

}  // namespace math
}  // namespace stan

#endif
