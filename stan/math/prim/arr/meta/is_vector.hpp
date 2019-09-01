#ifndef STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {

namespace internal {

/**
 * This underlying implimentation is used when the type is not an std vector.
 */
template <typename T>
struct is_std_vector_impl : std::false_type {};

/**
 * This specialization implimentation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_std_vector_impl<std::vector<Args...>> : std::true_type {};

}  // namespace internal

/**
 * Checks if the decayed type of T is a standard vector.
 */
template <typename T>
struct is_std_vector<
    T, std::enable_if_t<internal::is_std_vector_impl<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan
#endif
