#ifndef STAN_MATH_PRIM_FUN_SELECT_HPP
#define STAN_MATH_PRIM_FUN_SELECT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * If first argument is true return the second argument,
 * else return the third argument.
 *
 * `select(c, y1, y0) = c ? y1 : y0`.
 *
 * @tparam T_true A stan `Scalar` type
 * @tparam T_false A stan `Scalar` type
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_true, typename T_false,
          typename ReturnT = return_type_t<T_true, T_false>,
          require_all_stan_scalar_t<T_true, T_false>* = nullptr>
inline ReturnT select(const bool c, const T_true y_true,
                      const T_false y_false) {
  return c ? ReturnT(y_true) : ReturnT(y_false);
}

/**
 * If first argument is true return the second argument,
 * else return the third argument. Eigen expressions are
 * evaluated so that the return type is the same for both branches.
 *
 * Both containers must have the same plain type. The scalar type
 * of the return is determined by the return_type_t<> type trait.
 *
 * Overload for use with two containers.
 *
 * @tparam T_true A container of stan `Scalar` types
 * @tparam T_false A container of stan `Scalar` types
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <
    typename T_true, typename T_false,
    typename T_return = return_type_t<T_true, T_false>,
    typename T_true_plain = promote_scalar_t<T_return, plain_type_t<T_true>>,
    typename T_false_plain = promote_scalar_t<T_return, plain_type_t<T_false>>,
    require_all_container_t<T_true, T_false>* = nullptr,
    require_all_same_t<T_true_plain, T_false_plain>* = nullptr>
inline T_true_plain select(const bool c, T_true&& y_true, T_false&& y_false) {
  check_matching_dims("select", "left hand side", y_true, "right hand side",
                      y_false);
  return c ? T_true_plain(std::forward<T_true>(y_true))
           : T_true_plain(std::forward<T_false>(y_false));
}

/**
 * If first argument is true return the second argument,
 * else return the third argument.
 *
 * Overload for use when the 'true' return is a container and the 'false'
 * return is a scalar
 *
 * If chosen, the scalar is returned as a container of the same size and
 * plain type as the provided argument. Consequently, any Eigen expressions are
 * evaluated.
 *
 * @tparam T_true A container of stan `Scalar` types
 * @tparam T_false A stan `Scalar` type
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_true, typename T_false,
          typename ReturnT = promote_scalar_t<return_type_t<T_true, T_false>,
                                              plain_type_t<T_true>>,
          require_container_t<T_true>* = nullptr,
          require_stan_scalar_t<T_false>* = nullptr>
inline ReturnT select(const bool c, const T_true& y_true,
                      const T_false& y_false) {
  if (c) {
    return y_true;
  } else {
    return apply_scalar_binary(
        y_true, y_false,
        [](const auto& y_true_inner, const auto& y_false_inner) {
          return y_false_inner;
        });
  }
}

/**
 * If first argument is true return the second argument,
 * else return the third argument.
 *
 * Overload for use when the 'true' return is a scalar and the 'false'
 * return is a container
 *
 * If chosen, the scalar is returned as a container of the same size and
 * plain type as the provided argument. Consequently, any Eigen expressions are
 * evaluated.
 *
 * @tparam T_true A stan `Scalar` type
 * @tparam T_false A container of stan `Scalar` types
 * @param c Boolean condition value.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_true, typename T_false,
          typename ReturnT = promote_scalar_t<return_type_t<T_true, T_false>,
                                              plain_type_t<T_false>>,
          require_stan_scalar_t<T_true>* = nullptr,
          require_container_t<T_false>* = nullptr>
inline ReturnT select(const bool c, const T_true y_true,
                      const T_false y_false) {
  if (c) {
    return apply_scalar_binary(
        y_true, y_false,
        [](const auto& y_true_inner, const auto& y_false_inner) {
          return y_true_inner;
        });
  } else {
    return y_false;
  }
}

/**
 * If first argument is true return the second argument,
 * else return the third argument. Overload for use with an Eigen
 * object of booleans, and two scalars.
 *
 * The chosen scalar is returned as an Eigen object of the same dimension
 * as the input Eigen argument
 *
 * @tparam T_bool type of Eigen boolean object
 * @tparam T_true A stan `Scalar` type
 * @tparam T_false A stan `Scalar` type
 * @param c Eigen object of boolean condition values.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_bool, typename T_true, typename T_false,
          require_eigen_array_vt<std::is_integral, T_bool>* = nullptr,
          require_all_stan_scalar_t<T_true, T_false>* = nullptr>
inline auto select(const T_bool c, const T_true y_true, const T_false y_false) {
  using ret_t = return_type_t<T_true, T_false>;
  return c
      .unaryExpr([y_true, y_false](bool cond) {
        return cond ? ret_t(y_true) : ret_t(y_false);
      })
      .eval();
}

/**
 * If first argument is true return the second argument,
 * else return the third argument. Overload for use with an Eigen
 * array of booleans, one Eigen array and a scalar as input.
 *
 * @tparam T_bool type of Eigen boolean object
 * @tparam T_true A stan `Scalar` type or Eigen Array type
 * @tparam T_false A stan `Scalar` type or Eigen Array type
 * @param c Eigen object of boolean condition values.
 * @param y_true Value to return if condition is true.
 * @param y_false Value to return if condition is false.
 */
template <typename T_bool, typename T_true, typename T_false,
          require_eigen_array_t<T_bool>* = nullptr,
          require_any_eigen_array_t<T_true, T_false>* = nullptr>
inline auto select(const T_bool c, const T_true y_true, const T_false y_false) {
  check_consistent_sizes("select", "boolean", c, "left hand side", y_true,
                         "right hand side", y_false);
  using ret_t = return_type_t<T_true, T_false>;
  return c.select(y_true, y_false).template cast<ret_t>().eval();
}

}  // namespace math
}  // namespace stan

#endif
