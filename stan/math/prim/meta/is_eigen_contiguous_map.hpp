#ifndef STAN_MATH_PRIM_META_IS_EIGEN_CONTIGUOUS_MAP_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_CONTIGUOUS_MAP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {


namespace internal {
template <typename T>
struct is_eigen_contiguous_map_impl
    : bool_constant<std::is_base_of<Eigen::MapBase<std::decay_t<T>>, std::decay_t<T>>::value> {};

}  // namespace internal

/**
 * Check if a type inherits from eigen @c MapBase and has inner and outer
 *   strides equal to zero.
 */
template <typename T, typename Enable = void>
struct is_eigen_contiguous_map
    : std::false_type {};

template <typename T>
struct is_eigen_contiguous_map<T,
 std::enable_if_t<internal::is_eigen_contiguous_map_impl<std::decay_t<T>::value> :
   bool_constant<T::OuterStrideAtCompileTime == 0 && T::InnerStrideAtCompileTime == 0> {};


}  // namespace stan

#endif
