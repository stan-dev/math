#ifndef STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <type_traits>

namespace stan {

template <typename T, typename = void>
struct is_eigen;

// Checks whether decayed type is inherits from EigenBase
template <typename T>
struct is_eigen<T, std::enable_if_t<std::is_base_of<Eigen::EigenBase<std::decay_t<T>>,
                std::decay_t<T>>::value>> : std::true_type {};

}  // namespace stan
#endif
