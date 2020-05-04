#ifndef STAN_MATH_PRIM_META_IS_EIGEN_CONTIGUOUS_MAP_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_CONTIGUOUS_MAP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

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
