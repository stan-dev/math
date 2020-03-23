#ifndef STAN_MATH_PRIM_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_VALUE_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * Primary template class for metaprogram to compute the type of
 * values stored in a container.
 *
 * tparam T type of container.
 */
template <typename T, typename = void>
struct value_type {
  using type = typename std::decay_t<T>;
};

/** \ingroup type_trait
 * Specialization for pointers returns the underlying value the pointer is
 * pointing to.
 */
template <typename T>
struct value_type<T, std::enable_if_t<std::is_pointer<T>::value>> {
  using type = typename value_type<std::remove_pointer<T>>::type;
};

/** \ingroup type_trait
 * Helper function for accessing underlying type.
 */
template <typename T>
using value_type_t = typename value_type<T>::type;

/** \ingroup macro_helpers
 * Adds unary require aliases that check the `value_type`.
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_UNARY_VALUE(check_type, checker, doxygen_group)       \
  /*! \ingroup doxygen_group */                                                \
  /*! \addtogroup check_type##_types */                                        \
  /*! @{ */                                                                    \
  /*! \brief Require value type satisfies checker */                           \
  template <typename T>                                                        \
  using require_##check_type##_vt                                              \
      = require_t<checker<value_type_t<std::decay_t<T>>>>;                     \
  /*! \brief Require value type does not satisfy checker */                    \
  template <typename T>                                                        \
  using require_not_##check_type##_vt                                          \
      = require_not_t<checker<value_type_t<std::decay_t<T>>>>;                 \
  /*! \brief Require all of the value types satisfy checker */                 \
  template <typename... Types>                                                 \
  using require_all_##check_type##_vt                                          \
      = require_all_t<checker<value_type_t<std::decay_t<Types>>>...>;          \
  /*! \brief Require any of the value types satisfy checker */                 \
  template <typename... Types>                                                 \
  using require_any_##check_type##_vt                                          \
      = require_any_t<checker<value_type_t<std::decay_t<Types>>>...>;          \
  /*! \brief Require none of the value types satisfy checker */                \
  template <typename... Types>                                                 \
  using require_all_not_##check_type##_vt                                      \
      = require_all_not_t<checker<value_type_t<std::decay_t<Types>>>...>;      \
  /*! \brief Require at least one of the value types do not satisfy checker */ \
  template <typename... Types>                                                 \
  using require_any_not_##check_type##_vt                                      \
      = require_any_not_t<checker<value_type_t<std::decay_t<Types>>>...>;      \
/*! @} */

/** \ingroup macro_helpers
 * Adds binary require aliases that check the `value_type`.
 * @param check_type The name of the type to check, used to define
 * `require_<check_type>_t`.
 * @param checker A struct that returns holds a boolean `value`
 * @param doxygen_group The doxygen group to add this requires to.
 */
#define STAN_ADD_REQUIRE_BINARY_VALUE(check_type, checker, doxygen_group)     \
  /*! \ingroup doxygen_group */                                               \
  /*! \addtogroup check_type##_types */                                       \
  /*! @{ */                                                                   \
  /*! \brief Require that value types of `T` and `S` satisfies checker */     \
  template <typename T, typename S>                                           \
  using require_##check_type##_vt = require_t<                                \
      checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>; \
  template <typename T,                                                       \
            typename S> /*! \brief Require value types of `T` and `S`         \
                           does not satisfy checker */                        \
  using require_not_##check_type##_vt                                         \
      = require_not_t<checker<value_type_t<std::decay_t<T>>,                  \
                              value_type_t<std::decay_t<S>>>>;                \
  template <typename T, typename... Types> /*! \brief Require value types of  \
                                              `T` and all of the `Types`      \
                                              satisfy checker */              \
  using require_all_##check_type##_vt                                         \
      = require_all_t<checker<value_type_t<std::decay_t<T>>,                  \
                              value_type_t<std::decay_t<Types>>>...>;         \
  template <typename T,                                                       \
            typename... Types> /*! \brief Require any of the value types of   \
                                  `Types` and `T` satisfy checker */          \
  using require_any_##check_type##_vt                                         \
      = require_any_t<checker<value_type_t<std::decay_t<T>>,                  \
                              value_type_t<std::decay_t<Types>>>...>;         \
  template <typename T,                                                       \
            typename... Types> /*! \brief None of the value types of `Types`  \
                                  and `T` satisfy checker */                  \
  using require_all_not_##check_type##_vt                                     \
      = require_all_not_t<checker<value_type_t<std::decay_t<T>>,              \
                                  value_type_t<std::decay_t<Types>>>...>;     \
  /*! \brief Any of the value types `Types` and `T` do not satisfy checker */ \
  template <typename T, typename... Types>                                    \
  using require_any_not_##check_type##_vt                                     \
      = require_any_not_t<checker<value_type_t<std::decay_t<T>>,              \
                                  value_type_t<std::decay_t<Types>>>...>;     \
  /*! @} */

}  // namespace stan

#endif
