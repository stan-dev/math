#ifndef STAN_MATH_PRIM_SCAL_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_EIGEN_HPP

#include <type_traits>

namespace stan {

// Checks whether decayed type is inherits from EigenBase
template <typename T, typename = void>
struct is_eigen : std::false_type {};

}  // namespace stan
#endif
