#ifndef STAN_MATH_PRIM_META_IS_FVAR_HPP
#define STAN_MATH_PRIM_META_IS_FVAR_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member function type which is defined to be false
 * as the primitive scalar types cannot be a stan::math::fvar type.
 */
template <typename T, typename = void>
struct is_fvar : std::false_type {};

/** \addtogroup require_stan_scalar
*  @{
*/
/**
 * Require that a type is `fvar`.
 */
STAN_ADD_REQUIRE_UNARY(fvar, is_fvar);
STAN_ADD_REQUIRE_UNARY_SCALAR(fvar, is_fvar);
STAN_ADD_REQUIRE_UNARY_VALUE(fvar, is_fvar);
/** @}*/
}  // namespace stan
#endif

return_type_t<T_x, T_k> f;
if (k < 50.0) {
  int p = value_of_rec(28.0 + 0.5 * k - 100.0 / (k + 5.0) + 1);
  f = von_mises_cdf_series(x, k, p);
  if (f < 0) {
    f = 0;
  }
  if (f > 1) {
    f = 1;
  }
} else {
  return von_mises_cdf_normalapprox(x, k);
}
return f;
}
