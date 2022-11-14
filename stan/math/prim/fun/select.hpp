#ifndef STAN_MATH_PRIM_FUN_SELECT_HPP
#define STAN_MATH_PRIM_FUN_SELECT_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the second argument if the first argument is true
 * and otherwise return the third argument.
 *
 * <code>select(c, y1, y0) = c ? y1 : y0</code>.
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

/**
 * Return the second argument if the first argument is true
 * and otherwise return the third argument. Overload for use with two Eigen
 * objects.
 *
 * @tparam T_true type of the true argument
 * @tparam T_false type of the false argument
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_true, typename T_false,
          typename T_return = return_type_t<T_true, T_false>,
          typename T_true_plain = promote_scalar_t<T_return, plain_type_t<T_true>>,
          typename T_false_plain = promote_scalar_t<T_return, plain_type_t<T_false>>,
          require_all_eigen_t<T_true, T_false>* = nullptr,
          require_all_same_t<T_true_plain, T_false_plain>* = nullptr>
inline T_true_plain select(const bool c, const T_true y_true, const T_false y_false) {
  if (c) {
    return y_true;
  } else {
    return y_false;
  }
}

/**
 * Return the second Eigen argument if the first argument is true
 * and otherwise return the second Eigen argument. Overload for use with one
 * scalar and one Eigen object. If chosen, the scalar is returned as an Eigen
 * object of the same size and type as the provided argument.
 *
 * @tparam T_true type of the true argument
 * @tparam T_false type of the false argument
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_true, typename T_false,
          typename ReturnT = promote_scalar_t<return_type_t<T_true, T_false>,
                                              plain_type_t<T_true>>,
          require_eigen_t<T_true>* = nullptr,
          require_stan_scalar_t<T_false>* = nullptr>
inline ReturnT select(const bool c, const T_true& y_true,
                      const T_false& y_false) {
  if (c) {
    return y_true;
  } else {
    return y_true.unaryExpr([&](auto&& y) { return y_false; });
  }
}

/**
 * Return the second Eigen argument if the first argument is true
 * and otherwise return the second Eigen argument. Overload for use with one
 * scalar and one Eigen object. If chosen, the scalar is returned as an Eigen
 * object of the same size and type as the provided argument.
 *
 * @tparam T_true type of the true argument
 * @tparam T_false type of the false argument
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_true, typename T_false,
          typename ReturnT = promote_scalar_t<return_type_t<T_true, T_false>,
                                              plain_type_t<T_false>>,
          require_stan_scalar_t<T_true>* = nullptr,
          require_eigen_t<T_false>* = nullptr>
inline ReturnT select(const bool c, const T_true y_true,
                      const T_false y_false) {
  if (c) {
    return y_false.unaryExpr([&](auto&& y) { return y_true; });
  } else {
    return y_false;
  }
}

/**
 * Return the second argument if the first argument is true
 * and otherwise return the third argument. Overload for use with an Eigen
 * object of booleans, and two scalars. The chosen scalar is returned as an
 * Eigen object of the same dimension as the input Eigen argument
 *
 * @tparam T_bool type of Eigen boolean object
 * @tparam T_true type of the true argument
 * @tparam T_false type of the false argument
 * @param c Eigen object of boolean condition values.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_bool, typename T_true, typename T_false,
          require_eigen_array_t<T_bool>* = nullptr,
          require_all_stan_scalar_t<T_true, T_false>* = nullptr>
inline auto select(const T_bool c, const T_true y_true, const T_false y_false) {
  return c.unaryExpr([&](bool cond) { return cond ? y_true : y_false; }).eval();
}

/**
 * Return the second argument if the first argument is true
 * and otherwise return the third argument. Overload for use with an Eigen
 * object of booleans, and at least one Eigen object as input.
 *
 * @tparam T_bool type of Eigen boolean object
 * @tparam T_true type of the true argument
 * @tparam T_false type of the false argument
 * @param c Eigen object of boolean condition values.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_bool, typename T_true, typename T_false,
          require_eigen_array_t<T_bool>* = nullptr,
          require_any_eigen_array_t<T_true, T_false>* = nullptr,
          require_any_stan_scalar_t<T_true, T_false>* = nullptr>
inline auto select(const T_bool c, const T_true y_true, const T_false y_false) {
  return c.select(y_true, y_false).eval();
}

}  // namespace math
}  // namespace stan

#endif
