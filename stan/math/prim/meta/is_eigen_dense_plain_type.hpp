#ifndef STAN_MATH_PRIM_META_IS_EIGEN_DENSE_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_DENSE_PLAIN_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/is_eigen_dense_base.hpp>
#include <stan/math/prim/meta/is_plain_type.hpp>
#include <type_traits>

namespace stan {

/**
 * Checks whether type T is derived from Eigen::DenseBase.
 * If true this will have a static member function named value with a type
 * of true, else value is false.
 * @tparam T Type to check if it is derived from `DenseBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_dense_plain_type
    : bool_constant<is_eigen_dense_base<T>::value && is_plain_type<T>::value> {
};

STAN_ADD_REQUIRE_UNARY(eigen_dense_plain_type, is_eigen_dense_plain_type,
                       require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_dense_plain_type, is_eigen_dense_plain_type,
                           require_eigens_types);

}  // namespace stan

#endif
