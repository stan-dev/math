#ifndef STAN_MATH_PRIM_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

/**
 * Check if type derives from EigenBase
 **/
template <typename T, typename Enable = void>
struct is_eigen : std::false_type {};

template <typename T>
struct is_eigen<T,
 std::enable_if_t<std::is_base_of<Eigen::EigenBase<typename std::decay_t<std::decay_t<T>>::PlainObject>, typename std::decay_t<std::decay_t<T>>::PlainObject>::value>> : std::true_type {};

template <typename T>
struct is_eigen<T,
std::enable_if_t<std::is_base_of<Eigen::EigenBase<typename std::decay_t<std::decay_t<T>>::MatrixType>, typename std::decay_t<std::decay_t<T>>::MatrixType>::value>> : std::true_type {};

template <typename T>
struct is_eigen<Eigen::EigenBase<T>, void>: std::true_type {};

}  // namespace stan
#endif
