#ifndef STAN_MATH_PRIM_META_IS_CONTAINER_OR_VAR_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_CONTAINER_OR_VAR_MATRIX_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/**
 * Deduces whether type is eigen matrix, standard vector, or var<Matrix>.
 * @tparam Container type to check
 */
template <typename Container>
using is_container_or_var_matrix
    = bool_constant<math::disjunction<is_container<Container>,
                                      is_var_matrix<Container>>::value>;

}  // namespace stan

#endif
