#ifndef STAN_MATH_PRIM_SCAL_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_SCALAR_TYPE_HPP

#include <type_traits>

namespace stan {
/**
 * Metaprogram structure to determine the base scalar type
 * of a template argument.
 *
 * <p>This base class should be specialized for structured types.
 *
 * @tparam T Type of object.
 */
template <typename T, typename = void>
struct scalar_type {
  typedef T type;
};

template <typename T>
struct scalar_type<T, std::enable_if_t<std::is_pointer<T>::value>> {
  typedef typename scalar_type<T>::type type;
};

template <typename T>
using scalar_type_t = typename scalar_type<T>::type;

}  // namespace stan
#endif
