#ifndef STAN_MATH_PRIM_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Base implementation to check whether a type is derived from EigenBase
 */
template <typename T, typename = void>
struct is_eigen : std::false_type {};

namespace internal {
/**
 * Underlying implementation to check if a type is derived from EigenBase
 */
template <typename T>
struct is_eigen_base {
  static std::false_type f(const void *);
  template <typename Derived>
  static std::true_type f(const Eigen::EigenBase<Derived> *);
  enum { value = decltype(f(std::declval<T *>()))::value };
};

}  // namespace internal

/**
 * Checks whether type T is derived from EigenBase. If true this will have a
 * static member function named value with a type of true, else value is false.
 */
template <typename T>
struct is_eigen<
    T, std::enable_if_t<internal::is_eigen_base<std::decay_t<T>>::value>>
    : std::true_type {};

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

STAN_ADD_REQUIRE_UNARY(eigen, is_eigen, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen, is_eigen, require_eigens_types);

namespace internal {
template <typename T>
struct is_eigen_matrix_impl : std::false_type {};
template <typename T, int R, int C>
struct is_eigen_matrix_impl<Eigen::Matrix<T, R, C>> : std::true_type {};
template <typename T>
struct is_eigen_matrix_impl<Eigen::SparseMatrix<T>> : std::true_type {};
}  // namespace internal

/**
 * Check if a type is an `Eigen::Matrix` or `Eigen::SparseMatrix`
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_matrix : internal::is_eigen_matrix_impl<std::decay_t<T>> {};

STAN_ADD_REQUIRE_UNARY(eigen_matrix, is_eigen_matrix, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_matrix, is_eigen_matrix, require_eigens_types);

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
