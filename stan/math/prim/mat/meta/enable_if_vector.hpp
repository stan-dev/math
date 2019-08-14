#ifndef STAN_MATH_PRIM_MAT_META_ENABLE_IF_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_META_ENABLE_IF_VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

template <typename T>
using enable_if_vector = std::enable_if_t<is_vector<T>::value>;

template <typename T>
using enable_if_not_vector = std::enable_if_t<!is_vector<T>::value>;

template <typename... Types>
using enable_if_all_vector
    = std::enable_if_t<math::conjunction<is_vector<Types>...>::value>;

template <typename... Types>
using enable_if_any_vector
    = std::enable_if_t<math::disjunction<is_vector<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_vector
    = std::enable_if_t<!math::conjunction<is_vector<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_vector
    = std::enable_if_t<!math::disjunction<is_vector<Types>...>::value>;

template <typename T>
using enable_if_std_vector = std::enable_if_t<is_std_vector<T>::value>;

template <typename T>
using enable_if_not_std_vector = std::enable_if_t<!is_std_vector<T>::value>;

template <typename... Types>
using enable_if_all_std_vector
    = std::enable_if_t<math::conjunction<is_std_vector<Types>...>::value>;

template <typename... Types>
using enable_if_any_std_vector
    = std::enable_if_t<math::disjunction<is_std_vector<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_std_vector
    = std::enable_if_t<!math::conjunction<is_std_vector<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_std_vector
    = std::enable_if_t<!math::disjunction<is_std_vector<Types>...>::value>;

}  // namespace stan
#endif
