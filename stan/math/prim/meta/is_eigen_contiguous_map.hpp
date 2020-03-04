#ifndef STAN_MATH_PRIM_META_IS_EIGEN_CONTIGUOUS_MAP_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_CONTIGUOUS_MAP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {


  namespace internal {
  template <typename T>
  struct is_eigen_contiguous_map_impl : std::false_type {};
  template <typename T, int Opts>
  struct is_eigen_contiguous_map_impl<Eigen::Map<T, Opts, Eigen::Stride<0, 0>>>
      : std::true_type {};

  }  // namespace internal

  template <typename T>
  struct is_eigen_contiguous_map
      : internal::is_eigen_contiguous_map_impl<std::decay_t<T>> {};


}  // namespace stan

#endif
