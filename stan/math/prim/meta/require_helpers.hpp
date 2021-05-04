#ifndef STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>

#include <type_traits>

namespace stan {

/**
 * If condition is true, template is enabled
 * @ingroup type_traits
 */
template <class Check>
using require_t = std::enable_if_t<Check::value>;

/**
 * If condition is false, template is disabled
 * @ingroup type_traits
 */
template <typename Check>
using require_not_t = std::enable_if_t<!Check::value>;

/**
 * If all conditions are true, template is enabled
 * Returns a type void if all conditions are true and otherwise fails.
 * @ingroup type_traits
 */
template <class... Checks>
using require_all_t = std::enable_if_t<math::conjunction<Checks...>::value>;

/**
 * If any condition is true, template is enabled.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 * @ingroup type_traits
 */
template <class... Checks>
using require_any_t = std::enable_if_t<math::disjunction<Checks...>::value>;

/**
 * If all conditions are false, template is enabled.
 *
 * Returns a type void if all of the conditions are false.
 * @ingroup type_traits
 */
template <class... Checks>
using require_all_not_t
    = std::enable_if_t<!math::disjunction<Checks...>::value>;

/**
 * If any condition is false, template is enabled.
 *
 * Returns a type void if any of the conditions are false.
 * @ingroup type_traits
 */
template <class... Checks>
using require_any_not_t
    = std::enable_if_t<!math::conjunction<Checks...>::value>;

/** \ingroup macro_helpers
 * Adds Unary require aliases.
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_UNARY(check_type, checker, doxygen_group)       \
  /*! \ingroup doxygen_group */                                          \
  /*! \defgroup check_type##_types check_type  */                        \
  /*! \addtogroup check_type##_types */                                  \
  /*! @{ */                                                              \
  /*! \brief Require type satisfies checker */                           \
  template <typename T>                                                  \
  using require_##check_type##_t = require_t<checker<std::decay_t<T>>>;  \
                                                                         \
  /*! \brief Require type does not satisfy checker */                    \
  template <typename T>                                                  \
  using require_not_##check_type##_t                                     \
      = require_not_t<checker<std::decay_t<T>>>;                         \
                                                                         \
  /*! \brief Require all of the types satisfy checker */                 \
  template <typename... Types>                                           \
  using require_all_##check_type##_t                                     \
      = require_all_t<checker<std::decay_t<Types>>...>;                  \
                                                                         \
  /*! \brief Require any of the types satisfy checker */                 \
  template <typename... Types>                                           \
  using require_any_##check_type##_t                                     \
      = require_any_t<checker<std::decay_t<Types>>...>;                  \
                                                                         \
  /*! \brief Require none of the types satisfy checker */                \
  template <typename... Types>                                           \
  using require_all_not_##check_type##_t                                 \
      = require_all_not_t<checker<std::decay_t<Types>>...>;              \
                                                                         \
  /*! \brief Require at least one of the types do not satisfy checker */ \
  template <typename... Types>                                           \
  using require_any_not_##check_type##_t                                 \
      = require_any_not_t<checker<std::decay_t<Types>>...>;              \
/*! @} */

/** \ingroup macro_helpers
 * Adds unary require aliases that check the `value_type`.
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_UNARY_INNER(check_type, checker, doxygen_group)       \
  /*! \ingroup doxygen_group */                                                \
  /*! \addtogroup check_type##_types */                                        \
  /*! @{ */                                                                    \
  /*! \brief Require value type satisfies checker */                           \
  template <typename T>                                                        \
  using require_vt_##check_type                                                \
      = require_t<checker<value_type_t<std::decay_t<T>>>>;                     \
                                                                               \
  /*! \brief Require value type does not satisfy checker */                    \
  template <typename T>                                                        \
  using require_not_vt_##check_type                                            \
      = require_not_t<checker<value_type_t<std::decay_t<T>>>>;                 \
                                                                               \
  /*! \brief Require all of the value types satisfy checker */                 \
  template <typename... Types>                                                 \
  using require_all_vt_##check_type                                            \
      = require_all_t<checker<value_type_t<std::decay_t<Types>>>...>;          \
                                                                               \
  /*! \brief Require any of the value types satisfy checker */                 \
  template <typename... Types>                                                 \
  using require_any_vt_##check_type                                            \
      = require_any_t<checker<value_type_t<std::decay_t<Types>>>...>;          \
                                                                               \
  /*! \brief Require none of the value types satisfy checker */                \
  template <typename... Types>                                                 \
  using require_all_not_vt_##check_type                                        \
      = require_all_not_t<checker<value_type_t<std::decay_t<Types>>>...>;      \
                                                                               \
  /*! \brief Require at least one of the value types do not satisfy checker */ \
  template <typename... Types>                                                 \
  using require_any_not_vt_##check_type                                        \
      = require_any_not_t<checker<value_type_t<std::decay_t<Types>>>...>;      \
                                                                               \
  /*! \brief Require scalar type satisfies checker */                          \
  template <typename T>                                                        \
  using require_st_##check_type                                                \
      = require_t<checker<scalar_type_t<std::decay_t<T>>>>;                    \
                                                                               \
  /*! \brief Require scalar type does not satisfy checker */                   \
  template <typename T>                                                        \
  using require_not_st_##check_type                                            \
      = require_not_t<checker<scalar_type_t<std::decay_t<T>>>>;                \
                                                                               \
  /*! \brief Require all of the scalar types satisfy checker */                \
  template <typename... Types>                                                 \
  using require_all_st_##check_type                                            \
      = require_all_t<checker<scalar_type_t<std::decay_t<Types>>>...>;         \
                                                                               \
  /*! \brief Require any of the scalar types satisfy checker */                \
  template <typename... Types>                                                 \
  using require_any_st_##check_type                                            \
      = require_any_t<checker<scalar_type_t<std::decay_t<Types>>>...>;         \
                                                                               \
  /*! \brief Require none of the scalar types satisfy checker */               \
  template <typename... Types>                                                 \
  using require_all_not_st_##check_type                                        \
      = require_all_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>;     \
                                                                               \
  /*! \brief Any of the scalar types do not satisfy checker */                 \
  template <typename... Types>                                                 \
  using require_any_not_st_##check_type                                        \
      = require_any_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>;     \
/*! @} */

/** \ingroup macro_helpers
 * Adds binary require aliases.
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_BINARY(check_type, checker, doxygen_group)          \
  /*! \ingroup doxygen_group */                                              \
  /*! \defgroup check_type##_types check_type  */                            \
  /*! \addtogroup check_type##_types */                                      \
  /*! @{ */                                                                  \
  /*! \brief Require types `T` and `S` satisfies checker */                  \
  template <typename T, typename S>                                          \
  using require_##check_type##_t                                             \
      = require_t<checker<std::decay_t<T>, std::decay_t<S>>>;                \
                                                                             \
  /*! \brief Require types `T` and `S` does not satisfy checker */           \
  template <typename T, typename S>                                          \
  using require_not_##check_type##_t                                         \
      = require_not_t<checker<std::decay_t<T>, std::decay_t<S>>>;            \
                                                                             \
  /*! \brief Require `T` and all of the `Types` satisfy checker */           \
  template <typename T, typename... Types>                                   \
  using require_all_##check_type##_t                                         \
      = require_all_t<checker<std::decay_t<T>, std::decay_t<Types>>...>;     \
                                                                             \
  /*! \brief Require any of the `Types` and `T` satisfy checker */           \
  template <typename T, typename... Types>                                   \
  using require_any_##check_type##_t                                         \
      = require_any_t<checker<std::decay_t<T>, std::decay_t<Types>>...>;     \
                                                                             \
  /*! \brief Require none of the `Types` and `T` satisfy checker */          \
  template <typename T, typename... Types>                                   \
  using require_all_not_##check_type##_t                                     \
      = require_all_not_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \
                                                                             \
  /*! \brief Any one of the `Types` and `T` do not satisfy */                \
  template <typename T, typename... Types>                                   \
  using require_any_not_##check_type##_t                                     \
      = require_any_not_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \
/*! @} */

/** \ingroup macro_helpers
 * Adds binary require aliases that check the `scalar_type`.
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_BINARY_INNER(check_type, checker, doxygen_group)      \
  /*! \ingroup doxygen_group */                                                \
  /*! \addtogroup check_type##_types */                                        \
  /*! @{ */                                                                    \
  /*! \brief Require that value types of `T` and `S` satisfies checker */      \
  template <typename T, typename S>                                            \
  using require_st_##check_type                                                \
      = require_t<checker<scalar_type_t<std::decay_t<T>>,                      \
                          scalar_type_t<std::decay_t<S>>>>;                    \
                                                                               \
  /*! \brief Require scalar types of `T` and `S` does not satisfy checker */   \
  template <typename T, typename S>                                            \
  using require_not_st_##check_type                                            \
      = require_not_t<checker<scalar_type_t<std::decay_t<T>>,                  \
                              scalar_type_t<std::decay_t<S>>>>;                \
                                                                               \
  /*! \brief All scalar types of `T` and all of the `Types` satisfy checker */ \
  template <typename T, typename... Types>                                     \
  using require_all_st_##check_type                                            \
      = require_all_t<checker<scalar_type_t<std::decay_t<T>>,                  \
                              scalar_type_t<std::decay_t<Types>>>...>;         \
                                                                               \
  /*! \brief Any of the scalar types of `Types` and `T` satisfy checker */     \
  template <typename T, typename... Types>                                     \
  using require_any_st_##check_type                                            \
      = require_any_t<checker<scalar_type_t<std::decay_t<T>>,                  \
                              scalar_type_t<std::decay_t<Types>>>...>;         \
                                                                               \
  /*! \brief None of the scalar types of `Types` and `T` satisfy checker */    \
  template <typename T, typename... Types>                                     \
  using require_all_not_st_##check_type                                        \
      = require_all_not_t<checker<scalar_type_t<std::decay_t<T>>,              \
                                  scalar_type_t<std::decay_t<Types>>>...>;     \
                                                                               \
  /*! \brief Any of the scalar types `Types` and `T` do not satisfy checker */ \
  template <typename T, typename... Types>                                     \
  using require_any_not_st_##check_type                                        \
      = require_any_not_t<checker<scalar_type_t<std::decay_t<T>>,              \
                                  scalar_type_t<std::decay_t<Types>>>...>;     \
                                                                               \
  /*! \brief Value types of `T` and `S` satisfies checker */                   \
  template <typename T, typename S>                                            \
  using require_vt_##check_type = require_t<                                   \
      checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>;  \
                                                                               \
  /*! \brief Value types of `T` and `S` does not satisfy checker */            \
  template <typename T, typename S>                                            \
  using require_not_vt_##check_type = require_not_t<                           \
      checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>;  \
                                                                               \
  /*! \brief Value types of `T` and all of the `Types` satisfy checker */      \
  template <typename T, typename... Types>                                     \
  using require_all_vt_##check_type                                            \
      = require_all_t<checker<value_type_t<std::decay_t<T>>,                   \
                              value_type_t<std::decay_t<Types>>>...>;          \
                                                                               \
  /*! \brief Any of the value types of `Types` and `T` satisfy checker */      \
  template <typename T, typename... Types>                                     \
  using require_any_vt_##check_type                                            \
      = require_any_t<checker<value_type_t<std::decay_t<T>>,                   \
                              value_type_t<std::decay_t<Types>>>...>;          \
                                                                               \
  /*! \brief None of the value types of `Types` and `T` satisfy checker */     \
  template <typename T, typename... Types>                                     \
  using require_all_not_vt_##check_type                                        \
      = require_all_not_t<checker<value_type_t<std::decay_t<T>>,               \
                                  value_type_t<std::decay_t<Types>>>...>;      \
                                                                               \
  /*! \brief Any of the value types `Types` and `T` do not satisfy checker */  \
  template <typename T, typename... Types>                                     \
  using require_any_not_vt_##check_type                                        \
      = require_any_not_t<checker<value_type_t<std::decay_t<T>>,               \
                                  value_type_t<std::decay_t<Types>>>...>;      \
/*! @} */

/**
 * Used as the base for checking whether a type is a container with
 * an underlying scalar type
 *
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container satisfies a particular type
 * check.
 * @tparam ValueCheck Templated struct or alias that returns a containers
 * inner type.
 * @tparam CheckType Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container's underlying type satisfies a
 * particular type check.
 *
 */
template <template <class...> class ContainerCheck,
          template <class...> class ValueCheck,
          template <class...> class TypeCheck, class... Check>
using container_type_check_base
    = bool_constant<math::conjunction<ContainerCheck<std::decay_t<Check>>...,
                                      TypeCheck<ValueCheck<Check>>...>::value>;

/** \ingroup macro_helpers
 * Adds container require aliases that check the containers inner type.
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_CONTAINER(check_type, checker, doxygen_group)         \
  /*! \ingroup doxygen_group */                                                \
  /*! \defgroup check_type##_types check_type  */                              \
  /*! \addtogroup check_type##_types */                                        \
  /*! @{ */                                                                    \
  /*! \brief Require type satisfies checker */                                 \
  /*! and value type satisfies `TypeCheck` */                                  \
  /*! @tparam TypeCheck The type trait to check the value type against*/       \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_##check_type##_vt = require_t<                                 \
      container_type_check_base<checker, value_type_t, TypeCheck, Check...>>;  \
                                                                               \
  /*! \brief Require type does not satisfy checker or */                       \
  /*! value type does not satisfy `TypeCheck` */                               \
  /*! @tparam TypeCheck The type trait to check the value type against*/       \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_not_##check_type##_vt = require_not_t<                         \
      container_type_check_base<checker, value_type_t, TypeCheck, Check...>>;  \
                                                                               \
  /*! \brief Require any of the types satisfy checker */                       \
  /*! and any of the value types satisfy `TypeCheck` */                        \
  /*! @tparam TypeCheck The type trait to check the value type against*/       \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_##check_type##_vt = require_any_t<                         \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require at least one of the types does not satisfy checker */     \
  /*! and none of the value types satisfy `TypeCheck` */                       \
  /*! @tparam TypeCheck The type trait to check the value type against*/       \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_not_##check_type##_vt = require_any_not_t<                 \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require all of the types satisfy checker */                       \
  /*! and all of the value types satisfy `TypeCheck` */                        \
  /*! @tparam TypeCheck The type trait to check the value type against*/       \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_##check_type##_vt = require_all_t<                         \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require none of the types satisfy checker */                      \
  /*! and none of the value types satisfy `TypeCheck` */                       \
  /*! @tparam TypeCheck The type trait to check the value type against*/       \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_not_##check_type##_vt = require_all_not_t<                 \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require type satisfies checker */                                 \
  /*! and scalar type satisfies `TypeCheck` */                                 \
  /*! @tparam TypeCheck The type trait to check the scalar type against*/      \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_##check_type##_st = require_t<                                 \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check...>>; \
                                                                               \
  /*! \brief Require type does not satisfy checker */                          \
  /*! or scalar type does not satisfy `TypeCheck` */                           \
  /*! @tparam TypeCheck The type trait to check the scalar type against*/      \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_not_##check_type##_st = require_not_t<                         \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check...>>; \
                                                                               \
  /*! \brief Require any of the types satisfy checker */                       \
  /*! and any scalar type satisfies `TypeCheck` */                             \
  /*! @tparam TypeCheck The type trait to check the scalar type against*/      \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_##check_type##_st = require_any_t<                         \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
                                                                               \
  /*! \brief Require at least one of the types does not satisfy checker */     \
  /*! and any scalar type does not satisfy `TypeCheck` */                      \
  /*! @tparam TypeCheck The type trait to check the scalar type against*/      \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_not_##check_type##_st = require_any_not_t<                 \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
                                                                               \
  /*! \brief Require all of the types does not satisfy checker */              \
  /*! and all scalar types satisfy `TypeCheck` */                              \
  /*! @tparam TypeCheck The type trait to check the scalar type against*/      \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_##check_type##_st = require_all_t<                         \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
                                                                               \
  /*! \brief Require none of the types satisfy checker */                      \
  /*! and none of the scalar types satisfy `TypeCheck` */                      \
  /*! @tparam TypeCheck The type trait to check the scalar type against*/      \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_not_##check_type##_st = require_all_not_t<                 \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
  /*! @} */

}  // namespace stan

#endif
