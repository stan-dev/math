#ifndef STAN_MATH_PRIM_META_CONTAINS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_META_CONTAINS_STD_VECTOR_HPP

#include <stan/math/prim/meta/contains_std_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {
template <typename... Ts>
struct contains_std_vector : std::false_type {};
}  // namespace stan

namespace stan {

template <typename T, typename... Ts>
struct contains_std_vector<std::vector<T>, Ts...> : std::true_type {};

template <typename T, typename... Ts>
struct contains_std_vector<T, Ts...> : contains_std_vector<Ts...> {};

}  // namespace stan
#endif
