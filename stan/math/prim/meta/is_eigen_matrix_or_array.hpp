#ifndef STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_OR_ARRAY_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_MATRIX_OR_ARRAY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/is_eigen_array.hpp>
#include <type_traits>

namespace stan {

/** \addtogroup type_trait
 *  @{
 */

/**
 * Check if a type satisfies either of @c is_eigen_matrix or @c is_eigen_array
 * @tparam T Type to check if it is derived from `EigenBase`
 */
template <typename T>
using is_eigen_matrix_or_array
    = math::disjunction<is_eigen_matrix<std::decay_t<T>>,
                        is_eigen_array<std::decay_t<T>>>;
/** @}*/
}  // namespace stan

#endif
