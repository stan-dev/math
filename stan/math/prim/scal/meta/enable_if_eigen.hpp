#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_EIGEN_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_EIGEN_HPP

#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_eigen = std::enable_if_t<is_eigen<T>::value>;

template <typename T>
using enable_if_not_eigen = std::enable_if_t<!is_eigen<T>::value>;

template <typename... Types>
using enable_if_all_eigen
    = std::enable_if_t<math::conjunction<is_eigen<Types>...>::value>;

template <typename... Types>
using enable_if_any_eigen
    = std::enable_if_t<math::disjunction<is_eigen<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_eigen
    = std::enable_if_t<!math::conjunction<is_eigen<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_eigen
    = std::enable_if_t<!math::disjunction<is_eigen<Types>...>::value>;

}  // namespace stan
#endif
