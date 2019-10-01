#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VECTOR_LIKE_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VECTOR_LIKE_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {

template <typename T>
using enable_if_vector_like = std::enable_if_t<is_vector_like<T>::value>;

template <typename T>
using enable_if_not_vector_like = std::enable_if_t<!is_vector_like<T>::value>;

template <typename... Types>
using enable_if_all_vector_like
    = std::enable_if_t<math::conjunction<is_vector_like<Types>...>::value>;

template <typename... Types>
using enable_if_any_vector_like
    = std::enable_if_t<math::disjunction<is_vector_like<Types>...>::value>;

template <typename... Types>
using enable_if_all_not_vector_like
    = std::enable_if_t<!math::conjunction<is_vector_like<Types>...>::value>;

template <typename... Types>
using enable_if_any_not_vector_like
    = std::enable_if_t<!math::disjunction<is_vector_like<Types>...>::value>;

}  // namespace stan
#endif
