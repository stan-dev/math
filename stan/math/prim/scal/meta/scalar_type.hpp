#ifndef STAN_MATH_PRIM_SCAL_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_SCALAR_TYPE_HPP

#include <stan/math/prim/scal/meta/value_type.hpp>
#include <type_traits>

namespace stan {

/**
 * End of recursion to determine base scalar type of a type
 * @tparam T the type.
 * @return The underlying base scalar type. If T is not a container then this
 * returns T
 */
template <typename T, typename = void>
struct scalar_type_base {
  using type = T;
};

/**
 * Recursivly goes through array attributes until base type is found
 * @tparam T the type to clean and return.
 */
template <typename T>
struct scalar_type_base<
    T, std::enable_if_t<std::is_const<T>::value || std::is_pointer<T>::value>> {
  using type = typename scalar_type_base<
      std::remove_pointer_t<std::remove_cv_t<T>>>::type;
};

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
  using type = typename scalar_type_base<T>::type;
};

template <typename T>
using scalar_type_t = typename scalar_type<T>::type;

template <typename T>
using scalar_type_decay_t = typename scalar_type<std::decay_t<T>>::type;

}  // namespace stan
#endif
