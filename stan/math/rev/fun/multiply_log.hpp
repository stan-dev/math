#ifndef STAN_MATH_REV_FUN_MULTIPLY_LOG_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_LOG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/elt_multiply.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the value of a*log(b).
 *
 * When both a and b are 0, the value returned is 0
 * and no gradients are accumulated.
 * The partial derivative with respect to a is log(b).
 * The partial derivative with respect to b is a/b.
 *
 * @param a First variable.
 * @param b Second variable.
 * @return Value of a*log(b)
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr,
          require_any_var_t<T1, T2>* = nullptr>
inline var multiply_log(const T1& a, const T2& b) {
  if (value_of(a) == 0.0 && value_of(b) == 0.0){
    return var(0.0);
  }
  return make_callback_var(multiply_log(value_of(a), value_of(b)),
                          [a, b](const auto& res) mutable {
    if constexpr (!is_constant<T1>::value && !is_constant<T2>::value) {
      a.adj() += res.adj() * log(b.val());
      b.adj() += res.adj() * a.val() / b.val();
    } else if constexpr (!is_constant<T1>::value) {
      a.adj() += res.adj() * log(b);
    } else {
      b.adj() += res.adj() * a / b.val();
    }
  });
}

/**
 * Return the elementwise product `a * log(b)`.
 *
 * Both `T1` and `T2` are matrices, and one of `T1` or `T2` must be a
 * `var_value`
 *
 * @tparam T1 Type of first argument
 * @tparam T2 Type of second argument
 * @param a First argument
 * @param b Second argument
 * @return elementwise product of `a` and `log(b)`
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
          require_any_rev_matrix_t<T1, T2>* = nullptr>
inline auto multiply_log(T1&& a, T2&& b) {
  check_matching_dims("multiply_log", "a", a, "b", b);
  arena_t<T1> arena_a = std::forward<T1>(a);
  arena_t<T2> arena_b = std::forward<T2>(b);
  using return_t
    = return_var_matrix_t<decltype(multiply_log(value_of(arena_a), value_of(arena_b))), T1, T2>;
  arena_t<return_t> res = multiply_log(value_of(arena_a), value_of(arena_b));

  if constexpr (is_not_constant_v<T1> && is_not_constant_v<T2>) {
    reverse_pass_callback(
        [res, arena_a, arena_b]() mutable {
          auto is_zero = (arena_a.val().array() == 0.0 && arena_b.val().array() == 0.0);
          arena_a.adj().array()
              += is_zero.select(0.0, res.adj().array() * arena_b.val().array().log());
          arena_b.adj().array() += is_zero.select(0.0, res.adj().array() * arena_a.val().array()
                                   / arena_b.val().array());
        });
  } else if constexpr (is_not_constant_v<T1>) {
        reverse_pass_callback(
        [res, arena_a, arena_b]() mutable {
          auto is_zero = (arena_a.val().array() == 0.0 && arena_b.array() == 0.0);
          arena_a.adj().array()
                += is_zero.select(0.0, res.adj().array() * arena_b.array().log());
                             });
  } else {
        reverse_pass_callback(
        [res, arena_a, arena_b]() mutable {
          auto is_zero = (arena_a.array() == 0.0 && arena_b.val().array() == 0.0);
          arena_b.adj().array() += is_zero.select(0.0, res.adj().array()
                                  * arena_a.array()
                                  / arena_b.val().array());
                             });
  }
  return res;
}

/**
 * Return the product `a * log(b)`.
 * In the case where b is a scalar and and element of `a` and `b` are zero
 * the value returned is 0 and no gradients are accumulated.
 * For `b`'s adjoint, this function can still return NaN as the adjoint
 * of `b` if `a` is nonzero anywhere, but `b` is zero. Likewise,
 * `a`'s adjoint can have undefined values if `a` is nonzero but `b` is zero.
 *
 * @tparam T1 Type of matrix argument
 * @tparam T2 Type of scalar argument
 * @param a Matrix argument
 * @param b Scalar argument
 * @return Product of `a` and `log(b)`
 */
template <typename T1, typename T2, require_rev_matrix_t<T1>* = nullptr,
          require_stan_scalar_t<T2>* = nullptr>
inline auto multiply_log(T1&& a, T2&& b) {
  arena_t<T1> arena_a = a;
  using return_t
    = return_var_matrix_t<decltype(multiply_log(value_of(arena_a), value_of(b))), T1, T2>;
  arena_t<return_t> res = multiply_log(value_of(arena_a), value_of(b));
  if constexpr (is_not_constant_v<T1> && is_not_constant_v<T2>) {
    reverse_pass_callback(
        [res, arena_a, b]() mutable {
          auto is_zero = ((arena_a.val().array() == 0.0) + (b.val() == 0.0) > 1);
          arena_a.adj().array() += is_zero.select(0.0, res.adj().array() * log(b.val()));
            b.adj() += is_zero.select(0.0, (res.adj().array() * arena_a.val().array()) / b.val()).sum();
        });
  } else if constexpr (is_not_constant_v<T1>) {
    reverse_pass_callback(
        [res, arena_a, b]() mutable {
          auto is_zero = ((arena_a.val().array() == 0.0) + (b == 0.0) > 1);
          arena_a.adj().array() += is_zero.select(0.0, res.adj().array() * log(b));
        });
  } else {
    reverse_pass_callback(
        [res, arena_a, b]() mutable {
          auto is_zero = ((arena_a.array() == 0.0) + (b.val() == 0.0) > 1);
          b.adj() += is_zero.select(0.0, (res.adj().array() * arena_a.val().array()) / b.val()).sum();
        });
  }
  return res;
}

/**
 * Return the product `a * log(b)`.
 *
 * @tparam T1 Type of scalar argument
 * @tparam T2 Type of matrix argument
 * @param a Scalar argument
 * @param b Matrix argument
 * @return Product of `a` and `log(b)`
 */
template <typename T1, typename T2, require_stan_scalar_t<T1>* = nullptr,
          require_rev_matrix_t<T2>* = nullptr>
inline auto multiply_log(T1&& a, T2&& b) {
  arena_t<T2> arena_b = std::forward<T2>(b);
  using return_t
    = return_var_matrix_t<decltype(multiply_log(value_of(a), value_of(arena_b))), T1, T2>;
  arena_t<return_t> res = multiply_log(value_of(a), value_of(arena_b));
  if constexpr (is_not_constant_v<T1> && is_not_constant_v<T2>) {
    reverse_pass_callback(
        [res, a, arena_b]() mutable {
          auto is_zero = ((a.val() == 0.0) + (arena_b.val().array() == 0.0) > 1);
          a.adj() += is_zero.select(0.0, res.adj().array() * arena_b.val().array().log()).sum();
          arena_b.adj().array()
              += is_zero.select(0.0, a.val() * res.adj().array() / arena_b.val().array());
        });
  } else if constexpr (is_not_constant_v<T1>) {
    reverse_pass_callback(
        [res, a, arena_b]() mutable {
          auto is_zero = ((a.val() == 0.0) + (arena_b.array() == 0.0) > 1);
          a.adj() += is_zero.select(0.0, res.adj().array() * arena_b.array().log()).sum();
        });
  } else {
    reverse_pass_callback(
        [res, a, arena_b]() mutable {
          auto is_zero = ((a == 0.0) + (arena_b.val().array() == 0.0) > 1);
          arena_b.adj().array()
              += is_zero.select(0.0, a * res.adj().array() / arena_b.val().array());
        });
  }
  return res;
}

}  // namespace math
}  // namespace stan
#endif
