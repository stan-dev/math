#ifndef STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_MAT_META_IS_EIGEN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_eigen.hpp>
#include <type_traits>

namespace stan {

namespace internal {
/*
 * Underlying implimenation to check if a type is derived from EigenBase
 */
template <typename T>
struct is_eigen_base
    : std::integral_constant<bool,
                             std::is_base_of<Eigen::EigenBase<std::decay_t<T>>,
                                             std::decay_t<T>>::value> {};
}  // namespace internal

/*
 * Checks whether type T is derived from EigenBase. If true this will have a static member function named value with a type of true, else value is false.
 */
template <typename T>
struct is_eigen<T, std::enable_if_t<internal::is_eigen_base<T>::value>>
    : std::true_type {};

}  // namespace stan
#endif
