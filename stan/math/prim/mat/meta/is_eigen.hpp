#ifndef STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <type_traits>

namespace stan {

// Checks whether decayed type is inherits from EigenBase
template <typename T>
struct is_eigen_decay
    : std::integral_constant<bool,
                             std::is_base_of<Eigen::EigenBase<std::decay_t<T>>,
                                             std::decay_t<T>>::value> {};

}  // namespace stan
#endif
