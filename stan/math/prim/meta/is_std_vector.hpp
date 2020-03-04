#ifndef STAN_MATH_PRIM_META_IS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_META_IS_STD_VECTOR_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace internal {

/** \ingroup type_trait
 * This underlying implimentation is used when the type is not an std vector.
 */
template <typename T>
struct is_std_vector_impl : std::false_type {};

/** \ingroup type_trait
 * This specialization implimentation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_std_vector_impl<std::vector<Args...>> : std::true_type {};

}  // namespace internal

/** \ingroup type_trait
 * Base implimentation for checking if type is std vector
 */
template <typename T, typename = void>
struct is_std_vector : std::false_type {};

/** \ingroup type_trait
 * Checks if the decayed type of T is a standard vector.
 */
template <typename T>
struct is_std_vector<
    T, std::enable_if_t<internal::is_std_vector_impl<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan

#endif
