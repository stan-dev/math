#ifndef STAN_MATH_PRIM_META_PROMOTE_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_PROMOTE_SCALAR_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Template metaprogram to calculate a type for converting a
 * convertible type.  This is the base case.
 *
 * @tparam T result scalar type.
 * @tparam S input type
 */
template <typename T, typename S, typename Enable = void>
struct promote_scalar_type {
  /**
   * The promoted type.
   */
  using type = T;
};

/**
 * Template metaprogram to calculate a type for a container whose
 * underlying scalar is converted from the second template
 * parameter type to the first.
 *
 * @tparam T result scalar type.
 * @tparam S input type
 */
template <typename T, typename S>
struct promote_scalar_type<T, std::vector<S>> {
  /**
   * The promoted type.
   */
  using type = std::vector<typename promote_scalar_type<T, S>::type>;
};

/**
 * Specialization for `var_value` when the type to convert to is a `var_value`.
 * @tparam T a `var_value`.
 * @tparam S A `var_value` whose template type is derived from `EigenBase`.
 *  This specialization holds type `S`.
 */
template <typename T, typename S>
struct promote_scalar_type<
    T, S, require_all_t<is_var<T>, is_var<S>, is_eigen<value_type_t<S>>>> {
  /**
   * The promoted type.
   */
  using type = std::decay_t<S>;
};

/**
 * Specialization for `var_value` when the type to convert to is arithmetic.
 * @tparam T an arithmetic type.
 * @tparam S A `var_value` whose template type is derived from `EigenBase`.
 *  This specialization the promoted value type of `S`.
 */
template <typename T, typename S>
struct promote_scalar_type<T, S,
                           require_all_t<std::is_arithmetic<T>, is_var<S>,
                                         is_eigen<value_type_t<S>>>> {
  /**
   * The promoted type.
   */
  using type = typename promote_scalar_type<T, value_type_t<S>>::type;
};

/**
 * Template metaprogram to calculate a type for a matrix, vector, row vector or
 * Eigen::Array whose underlying scalar is converted from the second template
 * parameter type to the first.
 *
 * @tparam T result scalar type.
 * @tparam S input matrix type
 */
template <typename T, typename S>
struct promote_scalar_type<T, S, require_eigen_t<S>> {
  /**
   * The promoted type.
   */
  using type = typename std::conditional<
      std::is_same<typename Eigen::internal::traits<std::decay_t<S>>::XprKind,
                   Eigen::MatrixXpr>::value,
      Eigen::Matrix<typename promote_scalar_type<T, typename S::Scalar>::type,
                    S::RowsAtCompileTime, S::ColsAtCompileTime>,
      Eigen::Array<typename promote_scalar_type<T, typename S::Scalar>::type,
                   S::RowsAtCompileTime, S::ColsAtCompileTime>>::type;
};

template <typename T, typename S>
using promote_scalar_t =
    typename promote_scalar_type<std::decay_t<T>, std::decay_t<S>>::type;

}  // namespace math
}  // namespace stan
#endif
