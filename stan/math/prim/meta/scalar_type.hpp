#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
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

/** \ingroup type_trait
 * Specialization of scalar_type for vector to recursively return the inner
 * scalar type.
 *
 * @tparam T type of standard vector
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::value_type>;
};

/** \ingroup type_trait
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type of matrix
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::Scalar>;
};

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

#define STAN_ADD_REQUIRE_UNARY_SCALAR(check_type, checker) \
template <typename T> \
using require_##check_type##_st = require_t<checker<scalar_type_t<std::decay_t<T>>>>; \
template <typename T> \
using require_not_##check_type##_st = require_not_t<checker<scalar_type_t<std::decay_t<T>>>>;\
template <typename... Types> \
using require_all_##check_type##_st = \
 require_all_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \
template <typename... Types> \
using require_any_##check_type##_st = \
 require_any_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \
template <typename... Types> \
using require_all_not_##check_type##_st \
    = require_all_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \
template <typename... Types> \
using require_any_not_##check_type##_st \
    = require_any_not_t<checker<scalar_type_t<std::decay_t<Types>>>...>; \


    #define STAN_ADD_REQUIRE_BINARY_SCALAR(check_type, checker) \
    template <typename T, typename S> \
    using require_##check_type##_st = require_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<S>>>>; \
    template <typename T, typename S> \
    using require_not_##check_type##_st = require_not_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<S>>>>;\
    template <typename T, typename... Types> \
    using require_all_##check_type##_st = \
     require_all_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \
    template <typename T, typename... Types> \
    using require_any_##check_type##_st = \
     require_any_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \
    template <typename T, typename... Types> \
    using require_all_not_##check_type##_st \
        = require_all_not_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \
    template <typename T, typename... Types> \
    using require_any_not_##check_type##_st \
        = require_any_not_t<checker<scalar_type_t<std::decay_t<T>>, scalar_type_t<std::decay_t<Types>>>...>; \

        
}  // namespace stan
#endif
