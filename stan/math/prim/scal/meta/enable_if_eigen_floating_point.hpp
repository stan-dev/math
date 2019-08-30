#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_EIGEN_FLOATING_POINT_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_EIGEN_FLOATING_POINT_HPP

#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
struct is_eigen_floating_point
    : std::integral_constant<
          bool, is_eigen<T>::value
                    && std::is_floating_point<scalar_type_decay_t<T>>::value> {
};

template <typename T>
using enable_if_eigen_floating_point
    = std::enable_if_t<is_eigen_floating_point<T>::value>;

template <typename T>
using enable_if_not_eigen_floating_point
    = std::enable_if_t<!is_eigen_floating_point<T>::value>;

template <typename... Types>
using enable_if_all_eigen_floating_point = std::enable_if_t<
    math::conjunction<is_eigen_floating_point<Types>...>::value>;

template <typename... Types>
using enable_if_any_eigen_floating_point = std::enable_if_t<
    math::disjunction<is_eigen_floating_point<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_eigen_floating_point = std::enable_if_t<
    !math::conjunction<is_eigen_floating_point<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_eigen_floating_point = std::enable_if_t<
    !math::disjunction<is_eigen_floating_point<Types>...>::value>;

template <typename T>
using eigen_floating_point_type
    = enable_if_eigen_floating_point<std::decay_t<T>>;

template <typename T>
using not_eigen_floating_point_type
    = enable_if_not_eigen_floating_point<std::decay_t<T>>;

template <typename... Types>
using all_eigen_floating_point_type
    = enable_if_all_eigen_floating_point<std::decay_t<Types>...>;

template <typename... Types>
using any_eigen_floating_point_type
    = enable_if_any_eigen_floating_point<std::decay_t<Types>...>;

template <typename... Types>
using not_all_eigen_floating_point_type
    = enable_if_all_not_eigen_floating_point<std::decay_t<Types>...>;

template <typename... Types>
using not_any_eigen_floating_point_type
    = enable_if_any_not_eigen_floating_point<std::decay_t<Types>...>;

}  // namespace stan
#endif
