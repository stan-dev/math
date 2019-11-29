#ifndef STAN_MATH_PRIM_SCAL_META_IS_COMPLEX_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_COMPLEX_HPP

#include <type_traits>
#include <complex>

namespace stan {

template <typename T, typename = void>
struct is_complex : std::false_type {};


namespace internal {

  /** \ingroup type_trait
   * This underlying implimentation is used when the type is not std::complex.
   */
  template <typename T>
  struct is_complex_impl : std::false_type {};

  /** \ingroup type_trait
   * This specialization implimentation has a static member named value equal
   * to true when the template type is std::copmlex.
   */
  template <typename... Args>
  struct is_complex_impl<std::complex<Args...>> : std::true_type {};

  }  // namespace internal

/** \ingroup type_trait
 * Checks if the decayed type of T is std::complex.
 */
template <typename T>
struct is_complex<
    T, std::enable_if_t<internal::is_complex_impl<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan
#endif
