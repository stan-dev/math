#ifndef STAN_MATH_PRIM_ARR_META_CONTAINS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_ARR_META_CONTAINS_STD_VECTOR_HPP

#include <type_traits>
#include <vector>

namespace stan {

template <typename...>
struct contains_std_vector;

template <>
struct contains_std_vector<> : std::false_type {};

template <typename T, typename... Ts>
struct contains_std_vector<std::vector<T>, Ts...> : std::true_type {};

template <typename T, typename... Ts>
struct contains_std_vector<T, Ts...> : contains_std_vector<Ts...> {};

}  // namespace stan
#endif
