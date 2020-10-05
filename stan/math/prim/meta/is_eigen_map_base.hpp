#ifndef STAN_MATH_PRIM_META_IS_EIGEN_MAP_BASE_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_MAP_BASE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

namespace internal {
  template <typename T>
  using eigen_map_read_only = Eigen::MapBase<T, Eigen::ReadOnlyAccessors>;
  template <typename T>
  using eigen_map_write = Eigen::MapBase<T, Eigen::WriteAccessors>;
  template <typename T>
  using eigen_map_direct_read = Eigen::MapBase<T, Eigen::DirectAccessors>;
  template <typename T>
  using eigen_map_direct_write = Eigen::MapBase<T, Eigen::DirectWriteAccessors>;
}
/**
 * Checks whether type T is derived from Eigen::MapBase.
 * If true this will have a static member function named value with a type
 * of true, else value is false.
 * @tparam T Type to check if it is derived from `MapBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_map_base
    : bool_constant<is_base_pointer_convertible<internal::eigen_map_read_only, T>::value ||
     is_base_pointer_convertible<internal::eigen_map_write, T>::value ||
     is_base_pointer_convertible<internal::eigen_map_direct_read, T>::value ||
     is_base_pointer_convertible<internal::eigen_map_direct_write, T>::value> {
};

STAN_ADD_REQUIRE_UNARY(eigen_map_base, is_eigen_map_base,
                       require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(eigen_map_base, is_eigen_map_base,
                           require_eigens_types);

}  // namespace stan

#endif
