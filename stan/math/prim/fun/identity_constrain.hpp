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
 * transform to the input and increments the log probability
 * reference with the log absolute Jacobian determinant.
 *
 * This method is used for constrained types where x and the constrain differ
 *  in type.
 *
 * @tparam Scalar type of scalar
 * @tparam Types Types to check for promotion rules
 * @param[in] x scalar
 * @param[in] args values to check for promotion rules.
 * @return transformed input
 */
template <typename Scalar, typename... Types,
          require_all_stan_scalar_t<Scalar, Types...>* = nullptr>
inline auto identity_constrain(Scalar&& x, Types&&... args) {
  return return_type_t<Scalar, Types...>(x);
}

}  // namespace math
}  // namespace stan

#endif
