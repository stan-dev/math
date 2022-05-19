#ifndef STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>

#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * If condition is true, template is enabled
 * @ingroup type_trait
 */
template <class Check>
using require_t = std::enable_if_t<Check::value>;

/** \ingroup type_trait
 * If condition is false, template is disabled
 * @ingroup type_trait
 */
template <typename Check>
using require_not_t = std::enable_if_t<!Check::value>;

/** \ingroup type_trait
 * If all conditions are true, template is enabled
 * Returns a type void if all conditions are true and otherwise fails.
 * @ingroup type_trait
 */
template <class... Checks>
using require_all_t = std::enable_if_t<math::conjunction<Checks...>::value>;

/** \ingroup type_trait
 * If any condition is true, template is enabled.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 * @ingroup type_trait
 */
template <class... Checks>
using require_any_t = std::enable_if_t<math::disjunction<Checks...>::value>;

/** \ingroup type_trait
 * If all conditions are false, template is enabled.
 *
 * Returns a type void if all of the conditions are false.
 * @ingroup type_trait
 */
template <class... Checks>
using require_all_not_t
    = std::enable_if_t<!math::disjunction<Checks...>::value>;

/** \ingroup type_trait
 * If any condition is false, template is enabled.
 *
 * Returns a type void if any of the conditions are false.
 * @ingroup type_trait
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
#define STAN_ADD_REQUIRE_UNARY(check_type, checker, doxygen_group)      \
  /*! \ingroup doxygen_group */                                         \
  /*! \defgroup check_type##_types check_type  */                       \
  /*! \addtogroup check_type##_types */                                 \
  /*! @{ */                                                             \
  /*! \brief checker##<T> is true */                                    \
  template <typename T>                                                 \
  using require_##check_type##_t = require_t<checker<std::decay_t<T>>>; \
                                                                        \
  /*! \brief `##checker##<T>` is false */                               \
  template <typename T>                                                 \
  using require_not_##check_type##_t                                    \
      = require_not_t<checker<std::decay_t<T>>>;                        \
                                                                        \
  /*! \brief conjunction<##checker##<Types>...> is true */              \
  template <typename... Types>                                          \
  using require_all_##check_type##_t                                    \
      = require_all_t<checker<std::decay_t<Types>>...>;                 \
                                                                        \
  /*! \brief disjunction<##checker##<Types>...> is true */              \
  template <typename... Types>                                          \
  using require_any_##check_type##_t                                    \
      = require_any_t<checker<std::decay_t<Types>>...>;                 \
                                                                        \
  /*! \brief `conjunction<##checker##<Types>...>` is false */           \
  template <typename... Types>                                          \
  using require_all_not_##check_type##_t                                \
      = require_all_not_t<checker<std::decay_t<Types>>...>;             \
                                                                        \
  /*! \brief disjunction<##checker##<Types>...> is false */             \
  template <typename... Types>                                          \
  using require_any_not_##check_type##_t                                \
      = require_any_not_t<checker<std::decay_t<Types>>...>;             \
/*! @} */

/** \ingroup macro_helpers
 * Adds unary require aliases that check the \ref stan::value_type .
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_UNARY_INNER(check_type, checker, doxygen_group)     \
  /*! \ingroup doxygen_group */                                              \
  /*! \addtogroup check_type##_types */                                      \
  /*! @{ */                                                                  \
  /*! \brief `##checker##<value_type_t<T>>` is true */                       \
  template <typename T>                                                      \
  using require_vt_##check_type                                              \
      = require_t<checker<value_type_t<std::decay_t<T>>>>;                   \
                                                                             \
  /*! \brief `##checker##<value_type_t<T>>` is false */                      \
  template <typename T>                                                      \
  using require_not_vt_##check_type                                          \
      = require_not_t<checker<value_type_t<std::decay_t<T>>>>;               \
                                                                             \
  /*! \brief `conjunction<##checker##<value_type_t<Types>>...>` is true */   \
  template <typename... Types>                                               \
  using require_all_vt_##check_type                                          \
      = require_all_t<checker<value_type_t<std::decay_t<Types>>>...>;        \
                                                                             \
  /*! \brief `disjunction<##checker##<value_type_t<Types>>...>` is true */   \
  template <typename... Types>                                               \
  using require_any_vt_##check_type                                          \
      = require_any_t<checker<value_type_t<std::decay_t<Types>>>...>;        \
                                                                             \
  /*! \brief `conjunction<##checker##<value_type_t<Types>>...>` is false */  \
  template <typename... Types>                                               \
  using require_all_not_vt_##check_type                                      \
      = require_all_not_t<checker<value_type_t<std::decay_t<Types>>>...>;    \
                                                                             \
  /*! \brief `disjunction<##checker##<value_type_t<Types>>...>` is false */  \
  template <typename... Types>                                               \
  using require_any_not_vt_##check_type                                      \
      = require_any_not_t<checker<value_type_t<std::decay_t<Types>>>...>;    \
                                                                             \
  /*! \brief `##checker##<scalar_type_t<T>>` is true */                      \
  template <typename T>                                                      \
  using require_st_##check_type                                              \
      = require_t<checker<scalar_type_t<std::decay_t<T>>>>;                  \
                                                                             \
  /*! \brief `##checker##<scalar_type_t<T>>` is false */                     \
  template <typename T>                                                      \
  using require_not_st_##check_type                                          \
      = require_not_t<checker<scalar_type_t<std::decay_t<T>>>>;              \
                                                                             \
  /*! \brief `conjunction<##checker##<scalar_type_t<Types>>...>` is true */  \
  template <typename... Types>                                               \
  using require_all_st_##check_type                                          \
      = require_all_t<checker<scalar_type_t<std::decay_t<Types>>>...>;       \
                                                                             \
  /*! \brief `disjunction<##checker##<scalar_type_t<Types>>...>` is true */  \
  template <typename... Types>                                               \
  using require_any_st_##check_type                                          \
      = require_any_t<checker<scalar_type_t<std::decay_t<Types>>>...>;       \
                                                                             \
  /*! \brief `conjunction<##checker##<scalar_type_t<Types>>...>` is false */ \
  template <typename... Types>                                               \
  using require_all_not_st_##check_type                                      \
      = require_all_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>;   \
                                                                             \
  /*! \brief `disjunction<##checker##<scalar_type_t<Types>>...>` is false */ \
  template <typename... Types>                                               \
  using require_any_not_st_##check_type                                      \
      = require_any_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>;   \
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
  /*! \brief `conjunction<##checker##<T>, ##checker##<S>> is true` */        \
  template <typename T, typename S>                                          \
  using require_##check_type##_t                                             \
      = require_t<checker<std::decay_t<T>, std::decay_t<S>>>;                \
                                                                             \
  /*! \brief `conjunction<##checker##<T>, ##checker##<S>> is true` */        \
  template <typename T, typename S>                                          \
  using require_not_##check_type##_t                                         \
      = require_not_t<checker<std::decay_t<T>, std::decay_t<S>>>;            \
                                                                             \
  /*! \brief `conjunction<##checker##<T>, ##checker##<S>> is true` */        \
  template <typename T, typename... Types>                                   \
  using require_all_##check_type##_t                                         \
      = require_all_t<checker<std::decay_t<T>, std::decay_t<Types>>...>;     \
                                                                             \
  /*! \brief `conjunction<##checker##<T>, ##checker##<S>> is true` */        \
  template <typename T, typename... Types>                                   \
  using require_any_##check_type##_t                                         \
      = require_any_t<checker<std::decay_t<T>, std::decay_t<Types>>...>;     \
                                                                             \
  /*! \brief `conjunction<##checker##<T>, ##checker##<S>> is true` */        \
  template <typename T, typename... Types>                                   \
  using require_all_not_##check_type##_t                                     \
      = require_all_not_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \
                                                                             \
  /*! \brief `conjunction<##checker##<T>, ##checker##<S>> is true` */        \
  template <typename T, typename... Types>                                   \
  using require_any_not_##check_type##_t                                     \
      = require_any_not_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \
/*! @} */

/** \ingroup macro_helpers
 * Adds binary require aliases that check the \ref stan::scalar_type .
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
  /*! \brief Require `##checker##<Check>` is true and */                       \
  /*! \brief `TypeCheck<value_type_t<Check>>` is true */                       \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck The type trait to check the */                         \
  /*!  \ref stan::value_type against */                                        \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_##check_type##_vt = require_t<                                 \
      container_type_check_base<checker, value_type_t, TypeCheck, Check...>>;  \
                                                                               \
  /*! \brief Require `##checker##<Check>` is false or */                       \
  /*! \brief `TypeCheck<value_type_t<Check>>` is false */                      \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::value_type against*/    \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_not_##check_type##_vt = require_not_t<                         \
      container_type_check_base<checker, value_type_t, TypeCheck, Check...>>;  \
                                                                               \
  /*! \brief Require `disjunction<##checker##<Check>...>` is true and */       \
  /*! \brief `disjunction<TypeCheck<value_type_t<Check>>...>` is true */       \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::value_type against*/    \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_##check_type##_vt = require_any_t<                         \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require `disjunction<##checker##<Check>...>` is false or */       \
  /*! \brief `disjunction<TypeCheck<value_type_t<Check>>...>` is false */      \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::value_type against*/    \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_not_##check_type##_vt = require_any_not_t<                 \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require `conjunction<##checker##<Check>...>` is true and */       \
  /*! \brief `conjunction<TypeCheck<value_type_t<Check>>...>` is true */       \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::value_type against*/    \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_##check_type##_vt = require_all_t<                         \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require `conjunction<##checker##<Check>...>` is false and */      \
  /*! \brief `conjunction<TypeCheck<value_type_t<Check>>...>` is false */      \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::value_type against*/    \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_not_##check_type##_vt = require_all_not_t<                 \
      container_type_check_base<checker, value_type_t, TypeCheck, Check>...>;  \
                                                                               \
  /*! \brief Require `##checker##<Check>` is true and */                       \
  /*! \brief `TypeCheck<scalar_type_t<Check>>` is true */                      \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::scalar_type against*/   \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_##check_type##_st = require_t<                                 \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check...>>; \
                                                                               \
  /*! \brief Require `##checker##<Check>` is false or */                       \
  /*! \brief `TypeCheck<scalar_type_t<Check>>` is false */                     \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::scalar_type against*/   \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_not_##check_type##_st = require_not_t<                         \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check...>>; \
                                                                               \
  /*! \brief Require `disjunction<##checker##<Check>...>` is true and */       \
  /*! \brief `disjunction<TypeCheck<scalar_type_t<Check>...>>` is true */      \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::scalar_type against*/   \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_##check_type##_st = require_any_t<                         \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
                                                                               \
  /*! \brief Require `disjunction<##checker##<Check>...>` is false or */       \
  /*! \brief `disjunction<TypeCheck<scalar_type_t<Check>>...>` is false */     \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::scalar_type against*/   \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_any_not_##check_type##_st = require_any_not_t<                 \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
                                                                               \
  /*! \brief Require `conjunction<##checker##<Check>...>` is true and */       \
  /*! \brief `conjunction<TypeCheck<scalar_type_t<Check>>...>` is true */      \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::scalar_type against*/   \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_##check_type##_st = require_all_t<                         \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
                                                                               \
  /*! \brief Require `conjunction<##checker##<Check>...>` is false and */      \
  /*! \brief `conjunction<TypeCheck<scalar_type_t<Check>>...>` is false */     \
  /*! @tparam Check the type to check */                                       \
  /*! @tparam TypeCheck type trait to check \ref stan::scalar_type against*/   \
  template <template <class...> class TypeCheck, class... Check>               \
  using require_all_not_##check_type##_st = require_all_not_t<                 \
      container_type_check_base<checker, scalar_type_t, TypeCheck, Check>...>; \
  /*! @} */

}  // namespace stan

#endif
