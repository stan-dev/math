#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_vector = std::enable_if_t<is_vector<std::decay_t<T>>::value>;

template <typename... Types>
using enable_if_all_vector = std::enable_if_t<
    math::conjunction<is_vector<std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_any_vector = std::enable_if_t<
    math::disjunction<is_vector<std::decay_t<Types>>...>::value>;

template <typename T>
using enable_if_not_vector
    = std::enable_if_t<!is_vector<std::decay_t<T>>::value>;

template <typename... Types>
using enable_if_all_not_vector = std::enable_if_t<
    !math::conjunction<is_vector<std::decay_t<Types>>...>::value>;

template <typename... Types>
using enable_if_any_not_vector = std::enable_if_t<
    !math::disjunction<is_vector<std::decay_t<Types>>...>::value>;

template <typename T>
using vector_type = enable_if_vector<std::decay_t<T>>;

template <typename T>
using not_vector_type = enable_if_not_vector<std::decay_t<T>>;

template <typename... Types>
using all_vector_type = enable_if_all_vector<std::decay_t<Types>...>;

template <typename... Types>
using any_vector_type = enable_if_any_vector<std::decay_t<Types>...>;

template <typename... Types>
using not_all_vector_type = enable_if_all_not_vector<std::decay_t<Types>...>;

template <typename... Types>
using not_any_vector_type = enable_if_any_not_vector<std::decay_t<Types>...>;

}  // namespace stan
#endif
