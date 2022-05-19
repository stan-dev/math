#ifndef STAN_MATH_PRIM_META_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_META_PLAIN_TYPE_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_detected.hpp>
#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Determines plain (non expression) type associated with \c T. For non \c Eigen
 * types it is the decayed input type.
 * For Example, `plain_type<decltype(Eigen::MatrixXd() + Eigen::MatrixXd())>`
 * aka `plain_type<Eigen::CwiseBinaryOp<scalar_sum_op<double, double>, MatrixXd, MatrixXd>>`'s
 * `type` will be `Eigen::MatrixXd`. `plain_type<double>`'s `type` will be a double.
 * while `plain_type<std::vector<double>>`'s `type` will be a `double`
 * @tparam T type to determine plain type of
 */
template <typename T, typename Enable = void>
struct plain_type {
  using type = std::decay_t<T>;
};

/** \ingroup type_trait
 * Determine the non-expression type of an Eigen object. See \ref stan::plain_type for examples.
 */
template <typename T>
using plain_type_t = typename plain_type<T>::type;

/**
 * Determines return type of calling \c .eval() on Eigen expression.
 *
 * If input type \c T is a plain type (\c plain_type_t<T> equals \c
 * std::decay<T>), than member \c type is defined as <code> const
 * plain_type_t<T>& </code>. Otherwise member \c type is defined as \c
 * plain_type_t<T>.
 *
 * @tparam T type to determine eval return type of
 */
template <typename T>
struct eval_return_type {
  using T1 = plain_type_t<T>;
  using type = std::conditional_t<std::is_same<std::decay_t<T>, T1>::value,
                                  const T1&, T1>;
};

template <typename T>
using eval_return_type_t = typename eval_return_type<T>::type;

/** \ingroup type_trait
 * Determines plain (non expression) type associated with \c T. For \c Eigen
 * expression it is a type the expression can be evaluated into.
 * See the docs for \ref stan::plain_type for more information.
 * @tparam T type to determine plain type of
 */
template <typename T>
struct plain_type<T, require_eigen_t<T>> {
  using type = typename std::decay_t<T>::PlainObject;
};

}  // namespace stan

#endif  // STAN_MATH_PRIM_META_PLAIN_TYPE_HPP
