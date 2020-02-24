#ifndef STAN_MATH_PRIM_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {

/**
 * Check if type derives from EigenBase
 **/
template <typename T>
struct is_eigen : bool_constant<
                             std::is_base_of<Eigen::EigenBase<std::decay_t<T>>,
                              std::decay_t<T>>::value> {};

}  // namespace internal


}  // namespace stan
#endif
