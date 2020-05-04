#ifndef STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_OR_ARRAY_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_OR_ARRAY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/is_eigen_array.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

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


}  // namespace stan
#endif
