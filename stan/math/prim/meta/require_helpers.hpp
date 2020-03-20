#ifndef STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>



namespace stan {

  /**
   * If condition is true, template is enabled
   */
  template <class Check>
  using require_t = std::enable_if_t<Check::value>;

  /**
   * If condition is false, template is disabled
   */
  template <typename Check>
  using require_not_t = std::enable_if_t<!Check::value>;

  /**
   * If all conditions are true, template is enabled
   * Returns a type void if all conditions are true and otherwise fails.
   */
  template <class... Checks>
  using require_all_t = std::enable_if_t<math::conjunction<Checks...>::value>;

  /**
   * If any condition is true, template is enabled.
   *
   * Returns a type void if any of the conditions are true and otherwise fails.
   */
  template <class... Checks>
  using require_any_t = std::enable_if_t<math::disjunction<Checks...>::value>;

  /**
   * If all conditions are false, template is enabled.
   *
   * Returns a type void if all of the conditions are false.
   */
  template <class... Checks>
  using require_all_not_t
      = std::enable_if_t<!math::disjunction<Checks...>::value>;

  /**
   * If any condition is false, template is enabled.
   *
   * Returns a type void if any of the conditions are false.
   */
  template <class... Checks>
  using require_any_not_t
      = std::enable_if_t<!math::conjunction<Checks...>::value>;


#define STAN_ADD_REQUIRE_UNARY(check_type, checker, doxygen_group) \
/*! \ingroup doxygen_group */\
/*! \defgroup check_type##_types check_type  */\
template <typename T> /*! \brief Require type satisfies checker \ingroup check_type##_types */\
using require_##check_type##_t = require_t<checker<std::decay_t<T>>>; \
template <typename T> /*! \brief Require type does not satisfy checker \ingroup check_type##_types */\
using require_not_##check_type##_t = require_not_t<checker<std::decay_t<T>>>; \
template <typename... Types> /*! \brief Require all of the types satisfy checker \ingroup check_type##_types */\
using require_all_##check_type##_t = \
 require_all_t<checker<std::decay_t<Types>>...>; \
template <typename... Types> /*! \brief Require any of the types satisfy checker \ingroup check_type##_types */\
using require_any_##check_type##_t = \
 require_any_t<checker<std::decay_t<Types>>...>; \
template <typename... Types> /*! \brief Require none of the types satisfy checker \ingroup check_type##_types */\
using require_all_not_##check_type##_t \
    = require_all_not_t<checker<std::decay_t<Types>>...>; \
template <typename... Types> /*! \brief Require at least one of the types do not satisfy checker \ingroup check_type##_types */\
using require_any_not_##check_type##_t \
= require_any_not_t<checker<std::decay_t<Types>>...>;


#define STAN_ADD_REQUIRE_BINARY(check_type, checker) \
template <typename T, typename S> \
using require_##check_type##_t = require_t<checker<std::decay_t<T>, std::decay_t<S>>>; \
template <typename T, typename S> \
using require_not_##check_type##_t = require_not_t<checker<std::decay_t<T>, std::decay_t<S>>>;\
template <typename T, typename... Types> \
using require_all_##check_type##_t = \
 require_all_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \
template <typename T, typename... Types> \
using require_any_##check_type##_t = \
 require_any_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \
template <typename T, typename... Types> \
using require_all_not_##check_type##_t \
    = require_all_not_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \
template <typename T, typename... Types> \
using require_any_not_##check_type##_t \
    = require_any_not_t<checker<std::decay_t<T>, std::decay_t<Types>>...>; \


}  // namespace stan

#endif
