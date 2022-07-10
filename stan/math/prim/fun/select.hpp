#ifndef STAN_MATH_PRIM_FUN_SELECT_HPP
#define STAN_MATH_PRIM_FUN_SELECT_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the second argument if the first argument is true
 * and otherwise return the second argument.
 *
 * <p>This is just a convenience method to provide a function
 * with the same behavior as the built-in ternary operator.
 * In general, this function behaves as if defined by
 *
 * <p><code>select(c, y1, y0) = c ? y1 : y0</code>.
 *
 * @tparam T_true type of the true argument
 * @tparam T_false type of the false argument
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_true, typename T_false,
          require_all_stan_scalar_t<T_true, T_false>* = nullptr>
inline auto select(const bool c, const T_true y_true, const T_false y_false) {
  return c ? y_true : y_false;
}

template <typename T_true, typename T_false,
          require_all_eigen_t<T_true, T_false>* = nullptr>
inline auto select(const bool c, const T_true y_true, const T_false y_false) {
  using return_scalar_t = return_type_t<T_true, T_false>;
  if (c) {
    return y_true.template cast<return_scalar_t>().eval();
  }

  return y_false.template cast<return_scalar_t>().eval();
}

template <typename T_true, typename T_false,
          require_eigen_t<T_true>* = nullptr,
          require_stan_scalar_t<T_false>* = nullptr>
inline plain_type_t<T_true> select(const bool c, const T_true& y_true,
                                    const T_false& y_false) {
  if (c) {
    return y_true;
  }

  return y_true.unaryExpr([&](auto&& y){ return y_false; });
}


template <typename T_true, typename T_false,
          require_stan_scalar_t<T_true>* = nullptr,
          require_eigen_t<T_false>* = nullptr>
inline plain_type_t<T_false> select(const bool c, const T_true y_true, const T_false y_false) {
  if (c) {
  return y_false.unaryExpr([&](auto&& y){ return y_true; });
  }

  return y_false;
}

template <typename T_bool, typename T_true, typename T_false,
          require_eigen_array_t<T_bool>* = nullptr,
          require_all_stan_scalar_t<T_true, T_false>* = nullptr>
inline auto select(const T_bool c, const T_true y_true,
                    const T_false y_false) {
  return c.unaryExpr([&](bool cond){
    return cond ? y_true : y_false;
  }).eval();
}

template <typename T_bool, typename T_true, typename T_false,
          require_eigen_array_t<T_bool>* = nullptr,
          require_any_eigen_array_t<T_true, T_false>* = nullptr>
inline auto select(const T_bool c,
                                              const T_true y_true,
                                              const T_false y_false) {
  return c.select(y_true, y_false).eval();
}


}  // namespace math
}  // namespace stan

#endif
