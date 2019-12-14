#ifndef STAN_MATH_PRIM_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_META_SCALAR_TYPE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * End of recursion to determine base scalar type of a type
 * The underlying base scalar type. If T is not a container then this
 * has a static member nameed type with the type T.
 * @tparam T the type.
 */
template <typename T, typename = void>
struct scalar_type_base {
  using type = T;
};

/** \ingroup type_trait
 * Metaprogram structure to determine the base scalar type
 * of a template argument.
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T Type of object.
 */
template <typename T, typename = void>
struct scalar_type {
  using type = typename scalar_type_base<std::remove_cv_t<T>>::type;
};

template <typename T>
using scalar_type_t = typename scalar_type<T>::type;

template <typename T>
using scalar_type_decay_t = typename scalar_type<std::decay_t<T>>::type;

/** \ingroup type_trait
 * Specialization of scalar_type for vector to recursivly return the inner
 * scalar type.
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::value_type>;
};

/** \ingroup type_trait
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type of matrix.
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::Scalar>;
};
}  // namespace stan
#endif
