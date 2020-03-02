#ifndef STAN_MATH_PRIM_META_IS_COMPLEX_HPP
#define STAN_MATH_PRIM_META_IS_COMPLEX_HPP

#include <complex>
#include <type_traits>

namespace stan {
namespace internal {

/**
 * Provides a member constant `value` which is equal to `true` if
 * `T` is an instance of `std::complex` and `false` otherwise.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_complex_impl : std::false_type {};

template <typename... Ts>
struct is_complex_impl<std::complex<Ts...>> : std::true_type {};

}  // namespace internal

/**
 * If `T` is an arithmetic type (that is, an instance of
 * `std::complex`) or a cv-qualified version thereof, provides the
 * member constant `value` equal `true`; for any other type the value is
 * `false`.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<
    T, std::enable_if_t<internal::is_complex_impl<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan

#endif
