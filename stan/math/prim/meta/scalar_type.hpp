#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * Metaprogram structure to determine the base scalar type
 * of a template argument.
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T type of non-container
 */
template <typename T, typename = void>
struct scalar_type {
  using type = std::decay_t<T>;
};

template <typename T>
using scalar_type_t = typename scalar_type<T>::type;






/** \ingroup macro_helpers
 * Adds unary require aliases that check the `scalar_type`.
 * @param check_type The name of the type to check, used to define `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_UNARY_SCALAR(check_type, checker, doxygen_group) \
/*! \ingroup doxygen_group */\
/*! \addtogroup check_type##_types */ \
/*! @{ */ \
/*! \brief Require scalar type satisfies checker */\
template <typename T> \
using require_##check_type##_st = require_t<checker<scalar_type_t<std::decay_t<T>>>>; \
template <typename T> \
/*! \brief Require scalar type does not satisfy checker */\
using require_not_##check_type##_st = require_not_t<checker<scalar_type_t<std::decay_t<T>>>>;\
/*! \brief Require all of the scalar types satisfy checker */\
template <typename... Types> \
using require_all_##check_type##_st = \
 require_all_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \
/*! \brief Require any of the scalar types satisfy checker */\
template <typename... Types> \
using require_any_##check_type##_st = \
 require_any_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \
/*! \brief Require none of the scalar types satisfy checker */ \
template <typename... Types> \
using require_all_not_##check_type##_st \
    = require_all_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \
/*! \brief Require at least one of the scalar types do not satisfy checker */ \
template <typename... Types> \
using require_any_not_##check_type##_st \
    = require_any_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \
/*! @} */

/** \ingroup macro_helpers
 * Adds binary require aliases that check the `scalar_type`.
 * @param check_type The name of the type to check, used to define `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_BINARY_SCALAR(check_type, checker, doxygen_group) \
/*! \ingroup doxygen_group */\
/*! \addtogroup check_type##_types */ \
/*! @{ */ \
/*! \brief Require that value types of `T` and `S` satisfies checker */\
template <typename T, typename S> \
using require_##check_type##_st = require_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<S>>>>; \
template <typename T, typename S> \
/*! \brief Require scalar types of `T` and `S` does not satisfy checker */\
using require_not_##check_type##_st = require_not_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<S>>>>;\
template <typename T, typename... Types> \
/*! \brief All scalar types of `T` and all of the `Types` satisfy checker */\
using require_all_##check_type##_st = \
 require_all_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \
template <typename T, typename... Types> \
/*! \brief Any of the scalar types of `Types` and `T` satisfy checker */\
using require_any_##check_type##_st = \
 require_any_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \
template <typename T, typename... Types> \
/*! \brief None of the scalar types of `Types` and `T` satisfy checker */ \
using require_all_not_##check_type##_st \
    = require_all_not_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \
/*! \brief Any of the scalar types `Types` and `T` do not satisfy checker */ \
template <typename T, typename... Types> \
using require_any_not_##check_type##_st \
    = require_any_not_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \
/*! @} */



}  // namespace stan
#endif
