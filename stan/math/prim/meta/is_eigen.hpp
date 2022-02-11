#ifndef STAN_MATH_PRIM_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen_matrix_base.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {

/**
 * Check if type derives from `EigenBase`
 * @tparam T Type to check if it is derived from `EigenBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 **/
template <typename T>
struct is_eigen
    : bool_constant<is_base_pointer_convertible<Eigen::EigenBase, T>::value> {};

/**
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::Scalar>;
};

/**
 * Template metaprogram defining the type of values stored in an
 * Eigen matrix, vector, or row vector.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = typename std::decay_t<T>::Scalar;
};

STAN_ADD_REQUIRE_UNARY(eigen, is_eigen, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen, is_eigen, require_eigens_types);

/**
 * Check if a type is derived from `Eigen::ArrayBase`
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_array
    : bool_constant<is_base_pointer_convertible<Eigen::ArrayBase, T>::value> {};

STAN_ADD_REQUIRE_UNARY(eigen_array, is_eigen_array, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_array, is_eigen_array, require_eigens_types);

/**
 * Check if a type is derived from `Eigen::MatrixBase` or `Eigen::ArrayBase`
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
using is_eigen_matrix_or_array
    = math::disjunction<is_eigen_matrix_base<T>, is_eigen_array<T>>;

STAN_ADD_REQUIRE_UNARY(eigen_matrix_or_array, is_eigen_matrix_or_array,
                       require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_matrix_or_array, is_eigen_matrix_or_array,
                           require_eigens_types);

namespace internal {
template <typename T>
struct is_eigen_contiguous_map_impl : std::false_type {};
template <typename T, int Opts>
struct is_eigen_contiguous_map_impl<Eigen::Map<T, Opts, Eigen::Stride<0, 0>>>
    : std::true_type {};

}  // namespace internal

/**
 * Check if a type is an `Eigen::Map` with contiguous stride
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_contiguous_map
    : internal::is_eigen_contiguous_map_impl<std::decay_t<T>> {};

STAN_ADD_REQUIRE_UNARY(eigen_contiguous_map, is_eigen_contiguous_map,
                       require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_contiguous_map, is_eigen_contiguous_map,
                           require_eigens_types);

}  // namespace stan
#endif
