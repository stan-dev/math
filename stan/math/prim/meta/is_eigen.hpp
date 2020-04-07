#ifndef STAN_MATH_PRIM_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {

/** \addtogroup type_trait
 *  @{
 */

/**
 * Check if type derives from `EigenBase`
 * @tparam T Type to check if it is derived from `EigenBase`
 * @tparam Enable used for SFINAE deduction.
 **/
template <typename T, typename Enable = void>
struct is_eigen : std::false_type {};

template <typename T>
struct is_eigen<T, std::enable_if_t<std::is_base_of<
                       Eigen::EigenBase<typename std::decay_t<T>::PlainObject>,
                       typename std::decay_t<T>::PlainObject>::value>>
    : std::true_type {};

template <typename T>
struct is_eigen<T, std::enable_if_t<std::is_base_of<
                       Eigen::EigenBase<typename std::decay_t<T>::MatrixType>,
                       typename std::decay_t<T>::MatrixType>::value>>
    : std::true_type {};

template <typename T>
struct is_eigen<Eigen::EigenBase<T>, void> : std::true_type {};

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
 * Template metaprogram defining the type of values stored in an
 * Eigen matrix, vector, or row vector.
 *
 * @tparam T type of matrix.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = typename std::decay_t<T>::Scalar;
};
/** @}*/

STAN_ADD_REQUIRE_UNARY(eigen, is_eigen, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen, is_eigen, require_eigens_types);

namespace internal {
template <typename T>
struct is_eigen_array_impl : std::false_type {};
template <typename T, int R, int C>
struct is_eigen_array_impl<Eigen::Array<T, R, C>> : std::true_type {};
}  // namespace internal

/**
 * Check if a type is an `Eigen::Array`
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_array : internal::is_eigen_array_impl<std::decay_t<T>> {};

STAN_ADD_REQUIRE_UNARY(eigen_array, is_eigen_array, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_array, is_eigen_array, require_eigens_types);

/**
 * Check if a type is an `Eigen::Matrix` or `Eigen::SparseMatrix` or
 * `Eigen::Array`
 * @ingroup type_trait
 */
template <typename T>
using is_eigen_matrix_or_array
    = math::disjunction<is_eigen_matrix<T>, is_eigen_array<T>>;

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
