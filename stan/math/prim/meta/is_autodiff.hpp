#ifndef STAN_MATH_PRIM_META_IS_AUTODIFF_HPP
#define STAN_MATH_PRIM_META_IS_AUTODIFF_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <complex>
#include <type_traits>

namespace stan {

/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_autodiff_scalar
    : bool_constant<math::disjunction<is_var<std::decay_t<T>>,
                                      is_fvar<std::decay_t<T>>>::value> {};

template <typename T>
inline constexpr bool is_autodiff_scalar_v = is_autodiff_scalar<T>::value;

namespace internal {
/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct is_autodiff : is_autodiff_scalar<scalar_type_t<std::decay_t<T>>> {};

template <typename T>
struct is_autodiff<T, require_std_vector_t<T>>
    : bool_constant<is_autodiff<typename std::decay_t<T>::value_type>::value> {
};
template <typename T>
struct is_autodiff<T, require_eigen_t<T>>
    : bool_constant<is_autodiff<typename std::decay_t<T>::Scalar>::value> {};

}  // namespace internal

/**
 * Checks if decayed @ref scalar_type_t is a var or fvar
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct is_autodiff : internal::is_autodiff<T> {};

template <typename T>
inline constexpr bool is_autodiff_v = internal::is_autodiff<T>::value;

template <typename... Types>
inline constexpr bool is_all_autodiff_v = (is_autodiff_v<Types> && ...);

template <typename... Types>
inline constexpr bool is_any_autodiff_v = (is_autodiff_v<Types> || ...);

/*! \ingroup require_stan_scalar_real */
/*! \defgroup autodiff_types autodiff  */
/*! \addtogroup autodiff_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_autodiff_scalar */
/*! @tparam T the type to check */
template <typename T>
using require_autodiff_scalar_t
    = require_t<is_autodiff_scalar<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_autodiff_scalar */
/*! @tparam T the type to check */
template <typename T>
using require_not_autodiff_scalar_t
    = require_not_t<is_autodiff_scalar<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_autodiff_scalar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_autodiff_scalar_t
    = require_all_t<is_autodiff_scalar<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_autodiff_scalar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_autodiff_scalar_t
    = require_any_t<is_autodiff_scalar<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_autodiff_scalar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_autodiff_scalar_t
    = require_all_not_t<is_autodiff_scalar<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref
 * is_autodiff_scalar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_autodiff_scalar_t
    = require_any_not_t<is_autodiff_scalar<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \addtogroup autodiff_types */
/*! @{ */

/*! \brief Require value type does not satisfy @ref is_autodiff_scalar */
/*! @tparam T A type with a valid overload of @ref value_type available */
template <typename T>
using require_not_vt_autodiff_scalar
    = require_not_t<is_autodiff_scalar<value_type_t<std::decay_t<T>>>>;

/*! \brief Require none of the value types satisfy @ref is_autodiff_scalar */
/*! @tparam Types The types with a valid overload of @ref value_type available
 */
template <typename... Types>
using require_all_not_vt_autodiff_scalar = require_all_not_t<
    is_autodiff_scalar<value_type_t<std::decay_t<Types>>>...>;

/*! \brief Require scalar type satisfies @ref is_autodiff_scalar */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_st_autodiff_scalar
    = require_t<is_autodiff_scalar<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type does not satisfy @ref is_autodiff_scalar */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_not_st_autodiff_scalar
    = require_not_t<is_autodiff_scalar<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require any of the scalar types satisfy is_autodiff_scalar */
/*! @tparam Types The types with a valid overload of @ref scalar_type available
 */
template <typename... Types>
using require_any_st_autodiff_scalar
    = require_any_t<is_autodiff_scalar<scalar_type_t<std::decay_t<Types>>>...>;
/*! @} */

}  // namespace stan

#endif
