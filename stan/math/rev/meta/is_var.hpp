#ifndef STAN_MATH_REV_META_IS_VAR_HPP
#define STAN_MATH_REV_META_IS_VAR_HPP

#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/rev/core.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier is a var.
 */
template <typename T>
struct is_var<T,
              std::enable_if_t<std::is_same<math::var, std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan
#endif
