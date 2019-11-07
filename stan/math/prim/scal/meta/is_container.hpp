#ifndef STAN_MATH_PRIM_SCAL_META_IS_CONTAINER_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_CONTAINER_HPP

#include <stan/math/prim/scal/meta/bool_constant.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <type_traits>

namespace stan {

/**
 * Deduces whether type is eigen matrix or standard vector.
 * @tparam T type to check
 */
template <typename Container>
using is_container = bool_constant<
    math::disjunction<is_eigen<Container>, is_std_vector<Container>>::value>;
}  // namespace stan

#endif
