#ifndef STAN_MATH_PRIM_META_IS_COMPLEX_HPP
#define STAN_MATH_PRIM_META_IS_COMPLEX_HPP

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/is_autodiff.hpp>
#include <stan/math/prim/meta/base_type.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

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
 * If `T` is a complex type (that is, an instance of
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

/** \ingroup type_trait
 *
 * Template metaprogram defining the scalar type for values
 * stored in a complex number.
 *
 * @tparam T type of complex number
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_complex<T>::value>> {
  using type = std::complex<typename std::decay_t<T>::value_type>;
};

/**
 * Template metaprogram defining the base type for values
 * stored in a complex number.
 *
 * @tparam T type of complex number
 * @ingroup type_trait
 */
template <typename T>
struct base_type<T, std::enable_if_t<is_complex<T>::value>> {
  using type = base_type_t<typename std::decay_t<T>::value_type>;
};

/**
 * If `T` is a complex type with an inner arithmetic type (that is, an instance
 * of `std::complex<Arithmetic>`) or a cv-qualified version thereof, provides
 * the member constant `value` equal `true`; for any other type the value is
 * `false`.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct is_complex_arithmetic : std::false_type {};

template <typename T>
struct is_complex_arithmetic<
    T, std::enable_if_t<stan::math::conjunction<
           std::is_arithmetic<base_type_t<T>>,
           internal::is_complex_impl<std::decay_t<T>>>::value>>
    : std::true_type {};

/**
 * If `T` is a complex type with an inner autodiff type (that is, an instance of
 * `std::complex<var>` or `std::complex<fvar<T>>`) or a cv-qualified version
 * thereof, provides the member constant `value` equal `true`; for any other
 * type the value is `false`.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct is_complex_ad : std::false_type {};

template <typename T>
struct is_complex_ad<T, std::enable_if_t<stan::math::conjunction<
                            is_autodiff<base_type_t<T>>,
                            internal::is_complex_impl<std::decay_t<T>>>::value>>
    : std::true_type {};

/*! \ingroup require_stan_scalar_complex */
/*! \defgroup complex_types complex  */
/*! \addtogroup complex_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_complex */
/*! @tparam T the type to check */
template <typename T>
using require_complex_t = require_t<is_complex<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_complex */
/*! @tparam T the type to check */
template <typename T>
using require_not_complex_t = require_not_t<is_complex<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_complex */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_complex_t = require_all_t<is_complex<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_complex */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_complex_t = require_any_t<is_complex<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_complex */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_complex_t
    = require_all_not_t<is_complex<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_stan_scalar_complex */
/*! \addtogroup complex_types */
/*! @{ */

/*! \brief Require value type satisfies @ref is_complex */
/*! @tparam T A type with a valid overload of @ref value_type available */
template <typename T>
using require_vt_complex = require_t<is_complex<value_type_t<std::decay_t<T>>>>;

/*! \brief Require value type does not satisfy @ref is_complex */
/*! @tparam T A type with a valid overload of @ref value_type available */
template <typename T>
using require_not_vt_complex
    = require_not_t<is_complex<value_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type does not satisfy @ref is_complex */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_not_st_complex
    = require_not_t<is_complex<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require type satisfies @ref is_eigen */
/*! and scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class Check>
using require_complex_bt
    = require_all_t<is_complex<Check>, TypeCheck<base_type_t<Check>>>;

/*! @} */

/**
 * If the `value_type` of the type `T` is of type
 *  `std::complex` or a cv-qualified version thereof, provides the
 * member constant `value` equal `true`; for any other type the value is
 * `false`.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_vt_complex : is_complex<value_type_t<std::decay_t<T>>> {};

/**
 * If the `value_type` of the type `T` is not of type
 * `std::complex` or a cv-qualified version thereof, provides the
 * member constant `value` equal `true`; for any other type the value is
 * `false`.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_vt_not_complex
    : bool_constant<!is_complex<value_type_t<std::decay_t<T>>>::value> {};

}  // namespace stan

#endif
