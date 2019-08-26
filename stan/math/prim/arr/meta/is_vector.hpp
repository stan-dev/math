#ifndef STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {

namespace internal {
template <class T, typename = void>
struct is_same_vector : std::false_type {};

template <typename T>
using decay_vector = std::vector<typename std::decay_t<T>::value_type,
                                 typename std::decay_t<T>::allocator_type>;

template <typename T>
struct is_same_vector<
    T, std::enable_if_t<std::is_same<std::decay_t<T>, decay_vector<T>>::value>>
    : std::true_type {};

}  // namespace internal

template <typename T>
struct is_vector<
    T, std::enable_if_t<internal::is_same_vector<std::decay_t<T>>::value>>
    : std::true_type {};

template <typename T>
struct is_std_vector<
    T, std::enable_if_t<internal::is_same_vector<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan
#endif
