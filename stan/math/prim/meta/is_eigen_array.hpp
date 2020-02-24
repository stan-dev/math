#ifndef STAN_MATH_PRIM_META_IS_EIGEN_ARRAY_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_ARRAY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

/**
 * Check if a type derives from eigen ArrayBase
 */
template <typename T>
struct is_eigen_array : bool_constant<std::is_base_of<Eigen::ArrayBase<std::decay_t<T>>, std::decay_t<T>>::value> {};

}  // namespace stan

#endif
