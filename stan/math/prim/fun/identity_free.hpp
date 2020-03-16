#ifndef STAN_MATH_PRIM_FUN_IDENTITY_FREE_HPP
#define STAN_MATH_PRIM_FUN_IDENTITY_FREE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the result of applying the inverse of the identity
 * constraint transform to the input.
 *
 * <p>This function is a no-op and mainly useful as a placeholder
 * in auto-generated code.
 *
 * @tparam Scalar type of value
 * @param[in] x value
 * @return value
 */
template <typename Scalar, require_all_stan_scalar_t<Scalar>* = nullptr>
inline decltype(auto) identity_free(Scalar&& x) {
  return std::forward<Scalar>(x);
}

/**
 * Returns the result of applying the inverse of the identity
 * constraint transform to the input.
 *
 * This function is for N-ary op functions to make sure type promotion happens
 * correctly, thus allowing the use of auto.
 *
 *
 * @tparam Scalar type of value
 * @tparam Types The N-ary op argument types to deduce the return type from.
 * @param[in] x value
 * @return value promoted to the least upper bound of all input types.
 */
 template <typename Scalar, typename... Types,
           require_all_stan_scalar_t<Scalar, Types...>* = nullptr>
 inline auto identity_free(Scalar&& x, Types&&... args) {
   return return_type_t<Scalar, Types...>(x);
 }

 /**
  * Specialization for when a containers return type is the same as the input.
  *
  * This function is for N-ary op functions to make sure type promotion happens
  * correctly, thus allowing the use of auto.
  *
  * @tparam Container type of container
  * @tparam Types The N-ary op argument types to deduce the return type from.
  * @param[in] x value
  * @return value promoted to the least upper bound of all input types.
  */
 template <typename Container, typename... Types, require_container_t<Container>* = nullptr,
  require_same_t<return_type_t<Container, Types...>, return_type_t<Container>>* = nullptr>
 inline decltype(auto) identity_free(Container&& x, Types&&... args) {
   return std::forward<Container>(x);
 }

 /**
  * Returns the result of applying the inverse of the identity
  * constraint transform to the input.
  *
  * This function is for N-ary op functions to make sure type promotion happens
  * correctly, thus allowing the use of auto.
  *
  * @tparam EigT type derived from `Eigen::EigenBase`
  * @tparam Types The N-ary op argument types to deduce the return type from.
  * @param[in] x value
  * @return value promoted to the least upper bound of all input types.
  */
 template <typename EigT, typename... Types, require_eigen_t<EigT>* = nullptr,
 require_not_same_t<return_type_t<EigT, Types...>, return_type_t<EigT>>* = nullptr>
 inline auto identity_free(EigT&& x, Types&&... args) {
   return x.template cast<return_type_t<EigT, Types...>>().eval();
 }

 /**
  * Returns the result of applying the inverse of the identity
  * constraint transform to the input.
  *
  * This function is for N-ary op functions to make sure type promotion happens
  * correctly, thus allowing the use of auto.
  *
  *
  * @tparam StdVec a stanard vector.
  * @tparam Types The N-ary op argument types to deduce the return type from.
  * @param[in] x value
  * @return value promoted to the least upper bound of all input types.
  */
 template <typename StdVec, typename... Types, require_std_vector_t<StdVec>* = nullptr,
 require_not_same_t<return_type_t<StdVec, Types...>, return_type_t<StdVec>>* = nullptr>
 inline auto identity_free(StdVec&& x, Types&&... args) {
   return std::vector<return_type_t<StdVec, Types...>>(x.begin(), x.end());
 }

}  // namespace math
}  // namespace stan

#endif
