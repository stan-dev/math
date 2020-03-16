#ifndef STAN_MATH_PRIM_FUN_IDENTITY_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_IDENTITY_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of applying the identity constraint
 * transform to the input.
 *
 * <p>This method is effectively a no-op and is mainly useful as a
 * placeholder in auto-generated code.
 *
 * @tparam Scalar Type of scalar.
 * @param[in] x free scalar
 * @return transformed input
 */
template <typename Scalar, require_all_stan_scalar_t<Scalar>* = nullptr>
inline decltype(auto) identity_constrain(Scalar&& x) {
  return std::forward<Scalar>(x);
}

/**
 * Returns the result of applying the identity constraint
 * transform to the input.
 *
 * This method is used for constrained types where x and the constrain differ
 *  in type.
 *
 * @tparam Scalar type of scalar
 * @tparam Types Types to check for promotion rules
 * @param[in] x scalar
 * @param[in] args values to check for promotion rules.
 * @return promoted input
 */
template <typename Scalar, typename... Types,
          require_all_stan_scalar_t<Scalar, Types...>* = nullptr>
inline auto identity_constrain(Scalar&& x, Types&&... args) {
  return return_type_t<Scalar, Types...>(x);
}

/**
 * Returns the result of applying the identity constraint
 * transform to the input.
 *
 * This method is used for constrained types where x and the constrain do not
 * differ in type.
 *
 * @tparam Container Either a type derived from `Eigen::EigenBase` or a standard
 * vector.
 * @tparam Types Types to check for promotion rules
 * @param[in] x a container
 * @param[in] args values to check for promotion rules.
 * @return returns the input `x`
 */
template <typename Container, typename... Types,
          require_container_t<Container>* = nullptr,
          require_same_t<return_type_t<Container, Types...>,
                         return_type_t<Container>>* = nullptr>
inline decltype(auto) identity_constrain(Container&& x, Types&&... args) {
  return std::forward<Container>(x);
}

/**
 * Returns the result of applying the identity constraint
 * transform to the input.
 *
 * This method is used for constrained types where x and the constrain differ
 *  in type.
 *
 * @tparam EigT type derived from `Eigen::EigenBase`
 * @tparam Types Types to check for promotion rules
 * @param[in] x Eigen object
 * @param[in] args values to check for promotion rules.
 * @return promoted input
 */
template <typename EigT, typename... Types, require_eigen_t<EigT>* = nullptr,
          require_not_same_t<return_type_t<EigT, Types...>,
                             return_type_t<EigT>>* = nullptr>
inline auto identity_constrain(EigT&& x, Types&&... args) {
  return x.template cast<return_type_t<EigT, Types...>>().eval();
}

/**
 * Returns the result of applying the identity constraint
 * transform to the input.
 *
 * This method is used for constrained types where x and the constrain differ
 *  in type.
 *
 * @tparam StdVec type of standard vector
 * @tparam Types Types to check for promotion rules
 * @param[in] x standard vector
 * @param[in] args values to check for promotion rules.
 * @return promoted input
 */
template <typename StdVec, typename... Types,
          require_std_vector_t<StdVec>* = nullptr,
          require_not_same_t<return_type_t<StdVec, Types...>,
                             return_type_t<StdVec>>* = nullptr>
inline auto identity_constrain(StdVec&& x, Types&&... args) {
  return std::vector<return_type_t<StdVec, Types...>>(x.begin(), x.end());
}
}  // namespace math
}  // namespace stan

#endif
