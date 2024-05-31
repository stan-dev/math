#ifndef STAN_MATH_PRIM_META_IS_STRING_CONVERTIBLE_HPP
#define STAN_MATH_PRIM_META_IS_STRING_CONVERTIBLE_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <string>

namespace stan {

/**
 * Deduces whether type is convertible to string
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
using is_string_convertible = std::is_convertible<T, std::string>;

/*! \ingroup require_std */
/*! \defgroup string_convertible_types string_convertible  */
/*! \addtogroup string_convertible_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_string_convertible */
/*! @tparam T the type to check */
template <typename T>
using require_string_convertible_t
    = require_t<is_string_convertible<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_string_convertible */
/*! @tparam T the type to check */
template <typename T>
using require_not_string_convertible_t
    = require_not_t<is_string_convertible<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_string_convertible */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_string_convertible_t
    = require_all_t<is_string_convertible<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_string_convertible */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_string_convertible_t
    = require_any_t<is_string_convertible<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_string_convertible */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_string_convertible_t
    = require_all_not_t<is_string_convertible<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy
 * @ref is_string_convertible */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_string_convertible_t
    = require_any_not_t<is_string_convertible<std::decay_t<Types>>...>;
/*! @} */

}  // namespace stan

#endif
