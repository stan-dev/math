#ifndef STAN_MATH_PRIM_SCAL_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_EIGEN_HPP

#include <type_traits>

namespace stan {

/**
 * Base implimentation to check whether a type is derived from EigenBase
 */
template <typename T, typename = void>
struct is_eigen : std::false_type {};

template<class T>
constexpr bool is_eigen_v = is_eigen<T>::value;

}  // namespace stan
#endif
