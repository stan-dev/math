#ifndef STAN_MATH_PRIM_META_IS_EIGEN_MAP_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_MAP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {

namespace internal {
  template <typename T>
  using eigen_map_read_only = Eigen::MapBase<T, Eigen::ReadOnlyAccessors>;
  template <typename T>
  using eigen_map_write = Eigen::MapBase<T, Eigen::WriteAccessors>;
  template <typename T>
  using eigen_map_direct_read_only = Eigen::MapBase<T, Eigen::DirectAccessors>;
  template <typename T>
  using eigen_map_direct_write_only = Eigen::MapBase<T, Eigen::DirectWriteAccessors>;
}
/**
 * Check if type derives from `EigenBase`
 * @tparam T Type to check if it is derived from `EigenBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 **/
template <typename T>
struct is_eigen_map
    : bool_constant<is_base_pointer_convertible<internal::eigen_map_read_only, T>::value ||
       is_base_pointer_convertible<internal::eigen_map_write, T>::value ||
       is_base_pointer_convertible<internal::eigen_map_direct_read_only, T>::value ||
       is_base_pointer_convertible<internal::eigen_map_direct_write_only, T>::value> {};

STAN_ADD_REQUIRE_UNARY(eigen_map, is_eigen_map, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_map, is_eigen_map, require_eigens_types);

}  // namespace stan
#endif
