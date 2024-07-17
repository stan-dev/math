#ifndef STAN_MATH_REV_FUN_FMA_HPP
#define STAN_MATH_REV_FUN_FMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fma.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

/**
 * The fused multiply-add function for three variables (C99).
 * This function returns the product of the first two arguments
 * plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
inline var fma(const var& x, const var& y, const var& z) {
  return make_callback_var(fma(x.val(), y.val(), z.val()), [x, y, z](auto& vi) {
    x.adj() += vi.adj() * y.val();
    y.adj() += vi.adj() * x.val();
    z.adj() += vi.adj();
  });
}

/**
 * The fused multiply-add function for two variables and a value
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$.
 *
 * @tparam Tc type of the summand
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tc, require_arithmetic_t<Tc>* = nullptr>
inline var fma(const var& x, const var& y, Tc&& z) {
  return make_callback_var(fma(x.val(), y.val(), z), [x, y](auto& vi) {
    x.adj() += vi.adj() * y.val();
    y.adj() += vi.adj() * x.val();
  });
}

/**
 * The fused multiply-add function for a variable, value, and
 * variable (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tb, require_arithmetic_t<Tb>* = nullptr>
inline var fma(const var& x, Tb&& y, const var& z) {
  return make_callback_var(fma(x.val(), y, z.val()), [x, y, z](auto& vi) {
    x.adj() += vi.adj() * y;
    z.adj() += vi.adj();
  });
}

/**
 * The fused multiply-add function for a variable and two values
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The double-based version
 * <code>::%fma(double, double, double)</code> is defined in
 * <code>&lt;cmath&gt;</code>.
 *
 * The derivative is
 *
 * \f$\frac{d}{d x} (x * y) + z = y\f$.
 *
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tb, typename Tc, require_all_arithmetic_t<Tb, Tc>* = nullptr>
inline var fma(const var& x, Tb&& y, Tc&& z) {
  return make_callback_var(fma(x.val(), y, z),
                           [x, y](auto& vi) { x.adj() += vi.adj() * y; });
}

/**
 * The fused multiply-add function for a value, variable, and
 * value (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{d}{d y} (x * y) + z = x\f$, and
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tc type of the summand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tc, require_all_arithmetic_t<Ta, Tc>* = nullptr>
inline var fma(Ta&& x, const var& y, Tc&& z) {
  return make_callback_var(fma(x, y.val(), z),
                           [x, y](auto& vi) { y.adj() += vi.adj() * x; });
}

/**
 * The fused multiply-add function for two values and a variable,
 * and value (C99).  This function returns the product of the
 * first two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tb, require_all_arithmetic_t<Ta, Tb>* = nullptr>
inline var fma(Ta&& x, Tb&& y, const var& z) {
  return make_callback_var(fma(x, y, z.val()),
                           [z](auto& vi) { z.adj() += vi.adj(); });
}

/**
 * The fused multiply-add function for a value and two variables
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, require_arithmetic_t<Ta>* = nullptr>
inline var fma(Ta&& x, const var& y, const var& z) {
  return make_callback_var(fma(x, y.val(), z.val()), [x, y, z](auto& vi) {
    y.adj() += vi.adj() * x;
    z.adj() += vi.adj();
  });
}

namespace internal {

template <bool DoSum, typename T>
inline auto conditional_sum(T&& x) {
  if constexpr (DoSum) {
    return x.sum();
  } else {
    return std::forward<T>(x);
  }
}

template <typename T1, typename T2, typename T3, typename T4>
inline auto fma_reverse_pass(T1& arena_x, T2& arena_y, T3& arena_z, T4& ret) {
  return [arena_x, arena_y, arena_z, ret]() mutable {
    auto&& x_arr = as_array_or_scalar(arena_x);
    auto&& y_arr = as_array_or_scalar(arena_y);
    auto&& z_arr = as_array_or_scalar(arena_z);
    if constexpr (!is_constant_v<T1>) {
      x_arr.adj() += conditional_sum<is_stan_scalar_v<T1>>(ret.adj().array() * value_of(y_arr));
    }
    if constexpr (!is_constant_v<T2>) {
      y_arr.adj() += conditional_sum<is_stan_scalar_v<T2>>(ret.adj().array() * value_of(x_arr));
    }
    if constexpr (!is_constant_v<T3>) {
      z_arr.adj() += conditional_sum<is_stan_scalar_v<T3>>(ret.adj().array());
    }
  };
}

}  // namespace internal

/**
 * The fused multiply-add function for three variables (C99).
 * This function returns the product of the first two arguments
 * plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename T1, typename T2, typename T3,
          require_any_matrix_t<T1, T2, T3>* = nullptr,
          require_var_t<return_type_t<T1, T2, T3>>* = nullptr>
inline auto fma(T1&& x, T2&& y, T3&& z) {
  arena_t<T1> arena_x = std::forward<T1>(x);
  arena_t<T2> arena_y = std::forward<T2>(y);
  arena_t<T3> arena_z = std::forward<T3>(z);
  if constexpr (is_matrix_v<T1, T2>) {
    check_matching_dims("fma", "x", arena_x, "y", arena_y);
  }
  if constexpr (is_matrix_v<T1, T3>) {
    check_matching_dims("fma", "x", arena_x, "z", arena_z);
  }
  if constexpr (is_matrix_v<T2, T3>) {
    check_matching_dims("fma", "y", arena_y, "z", arena_z);
  }
  using inner_ret_type
      = decltype(fma(value_of(arena_x), value_of(arena_y), value_of(arena_z)));
  using ret_type = return_var_matrix_t<inner_ret_type, T1, T2, T3>;
  arena_t<ret_type> ret
      = fma(value_of(arena_x), value_of(arena_y), value_of(arena_z));
  reverse_pass_callback(
      internal::fma_reverse_pass(arena_x, arena_y, arena_z, ret));
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
