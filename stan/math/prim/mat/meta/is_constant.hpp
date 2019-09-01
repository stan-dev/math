#ifndef STAN_MATH_PRIM_MAT_META_IS_CONSTANT_HPP
#define STAN_MATH_PRIM_MAT_META_IS_CONSTANT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/is_constant.hpp>
#include <type_traits>

namespace stan {
/**
 * Defines a public enum named value and sets it to true
 * if the type of the elements in the provided Eigen Matrix
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 *
 * @tparam T type of the Eigen Matrix
 */
template <typename T>
struct is_constant<T, std::enable_if_t<is_eigen<T>::value>>
    : bool_constant<is_constant<typename std::decay_t<T>::Scalar>::value> {};

}  // namespace stan
#endif
