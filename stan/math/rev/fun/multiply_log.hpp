#ifndef STAN_MATH_REV_FUN_MULTIPLY_LOG_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_LOG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/elt_multiply.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/select.hpp>
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
  if (value_of(a) == 0.0 && value_of(b) == 0.0) {
    return var(0.0);
  }
  return make_callback_var(
      multiply_log(value_of(a), value_of(b)), [a, b](const auto& res) mutable {
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

namespace internal {
template <bool Cond, typename T>
inline auto conditional_sum(T&& x) {
  if constexpr (Cond) {
    return sum(std::forward<T>(x));
  } else {
    return std::forward<T>(x);
  }
}
}  // namespace internal

/**
 * Return the elementwise product `a * log(b)`.
 *
 * For each element of `a` and `b`, when `a[i]` and `b[i]` are 0,
 *  the value and adjoint returned are zero.
 * @tparam T1 Either a scalar or a matrix
 * @tparam T2 Either a scalar or a matrix
 * @param a First argument
 * @param b Second argument
 * @return elementwise product of `a` and `log(b)`
 */
template <typename T1, typename T2, require_any_matrix_t<T1, T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline auto multiply_log(T1&& a, T2&& b) {
  constexpr bool is_a_scalar = !is_matrix_v<T1>;
  constexpr bool is_b_scalar = !is_matrix_v<T2>;
  if constexpr (!is_a_scalar && !is_b_scalar) {
    check_matching_dims("multiply_log", "a", a, "b", b);
  }
  arena_t<T1> arena_a = std::forward<T1>(a);
  arena_t<T2> arena_b = std::forward<T2>(b);
  using return_t = return_var_matrix_t<
      decltype(multiply_log(value_of(arena_a), value_of(arena_b))), T1, T2>;
  arena_t<return_t> res = multiply_log(value_of(arena_a), value_of(arena_b));
  using internal::conditional_sum;
  if constexpr (is_not_constant_v<T1> && is_not_constant_v<T2>) {
    reverse_pass_callback([res, arena_a, arena_b]() mutable {
      auto arena_a_arr = as_array_or_scalar(arena_a);
      auto arena_b_arr = as_array_or_scalar(arena_b);
      auto res_arr = as_array_or_scalar(res);
      auto is_zero
          = ((arena_a_arr.val() == 0.0 + arena_b_arr.val() == 0.0) > 1);
      arena_a_arr.adj() += conditional_sum<is_a_scalar>(
          select(is_zero, 0.0, res_arr.adj() * log(arena_b_arr.val())));
      arena_b_arr.adj() += conditional_sum<is_b_scalar>(select(
          is_zero, 0.0, res_arr.adj() * arena_a_arr.val() / arena_b_arr.val()));
    });
  } else if constexpr (is_not_constant_v<T1>) {
    reverse_pass_callback([res, arena_a, arena_b]() mutable {
      auto arena_a_arr = as_array_or_scalar(arena_a);
      auto arena_b_arr = as_array_or_scalar(arena_b);
      auto res_arr = as_array_or_scalar(res);
      auto is_zero = ((arena_a_arr.val() == 0.0 + arena_b_arr == 0.0) > 1);
      arena_a_arr.adj() += conditional_sum<is_a_scalar>(
          select(is_zero, 0.0, res_arr.adj() * log(arena_b_arr)));
    });
  } else {
    reverse_pass_callback([res, arena_a, arena_b]() mutable {
      auto arena_a_arr = as_array_or_scalar(arena_a);
      auto arena_b_arr = as_array_or_scalar(arena_b);
      auto res_arr = as_array_or_scalar(res);
      auto is_zero = ((arena_a_arr == 0.0 + arena_b_arr.val() == 0.0) > 1);
      arena_b_arr.adj() += conditional_sum<is_b_scalar>(select(
          is_zero, 0.0, res_arr.adj() * arena_a_arr / arena_b_arr.val()));
    });
  }
  return res;
}

}  // namespace math
}  // namespace stan
#endif
