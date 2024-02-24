#ifndef STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP
#define STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_detected.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

namespace internal {
/**
 * @brief Used to detect if object has operator[](int) defined
 */
template <typename T>
using operator_bracket_t = decltype(std::declval<T>()[int{}]);
}  // namespace internal

/** \ingroup type_trait
 * Template metaprogram indicates whether a type is vector_like.
 *
 * A type is vector_like if an instance can be accessed like a
 * vector, i.e. square brackets.
 *
 * Access is_vector_like::value for the result.
 *
 * Default behavior is to use the is_vector template metaprogram.
 *
 * @tparam T Type to test
 */
template <typename T>
struct is_vector_like
    : bool_constant<is_detected<T, internal::operator_bracket_t>::value> {};

/*! \ingroup require_std */
/*! \defgroup vector_like_types vector_like  */
/*! \addtogroup vector_like_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_vector_like */
/*! @tparam T the type to check */
template <typename T>
using require_vector_like_t = require_t<is_vector_like<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_vector_like */
/*! @tparam T the type to check */
template <typename T>
using require_not_vector_like_t
    = require_not_t<is_vector_like<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_vector_like */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_vector_like_t
    = require_all_t<is_vector_like<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_vector_like */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_vector_like_t
    = require_any_t<is_vector_like<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_vector_like */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_vector_like_t
    = require_all_not_t<is_vector_like<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_vector_like
 */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_vector_like_t
    = require_any_not_t<is_vector_like<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_std */
/*! \defgroup vector_like_types vector_like  */
/*! \addtogroup vector_like_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_vector_like */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_vector_like for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_vector_like_vt
    = require_t<container_type_check_base<is_vector_like, value_type_t,
                                          TypeCheck, Check...>>;
/*! @} */

}  // namespace stan
#endif
