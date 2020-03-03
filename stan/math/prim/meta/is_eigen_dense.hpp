#ifndef STAN_MATH_PRIM_META_IS_EIGEN_DENSE_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_DENSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

  /**
   * Check if type derives from DenseBase
   **/
  template <typename T, typename Enable = void>
  struct is_eigen_dense : std::false_type {};

  template <typename T>
  struct is_eigen_dense<T,
   std::enable_if_t<std::is_base_of<Eigen::DenseBase<typename std::decay_t<T>::PlainObject>, typename std::decay_t<T>::PlainObject>::value>> : std::true_type {};

  template <typename T>
  struct is_eigen_dense<T,
  std::enable_if_t<std::is_base_of<Eigen::DenseBase<typename std::decay_t<T>::MatrixType>, typename std::decay_t<T>::MatrixType>::value>> : std::true_type {};

}  // namespace stan

#endif
