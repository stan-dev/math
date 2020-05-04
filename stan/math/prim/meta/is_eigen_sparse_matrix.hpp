#ifndef STAN_MATH_PRIM_META_IS_EIGEN_SPARSE_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_SPARSE_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

/**
 * Check if a type is derived from `Eigen::ArrayBase`
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_sparse_matrix
    : bool_constant<is_base_pointer_convertible<Eigen::SparseMatrixBase, T>::value> {};

STAN_ADD_REQUIRE_UNARY(eigen_sparse_matrix, is_eigen_sparse_matrix, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_sparse_matrix, is_eigen_sparse_matrix, require_eigens_types);

}  // namespace stan
#endif
