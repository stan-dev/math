#ifndef STAN_MATH_PRIM_META_BASE_TYPE_HPP
#define STAN_MATH_PRIM_META_BASE_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/**
 * Metaprogram structure to determine the base base type of a template
 * argument.  Qualifiers `const` and `volatile` are removed from all
 * types as are references.
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T type of non-container
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct base_type {
  using type = std::decay_t<T>;
};

template <typename T>
using base_type_t = typename base_type<T>::type;

/**
 * Specialization of base_type for vector to recursively return the inner
 * base type.
 *
 * @tparam T type of standard vector
 * @ingroup type_trait
 */
template <typename T>
struct base_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = base_type_t<typename std::decay_t<T>::value_type>;
};

/**
 * Template metaprogram defining the base base type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type of matrix
 * @ingroup type_trait
 */
template <typename T>
struct base_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = base_type_t<typename std::decay_t<T>::Scalar>;
};

}  // namespace stan
#endif
