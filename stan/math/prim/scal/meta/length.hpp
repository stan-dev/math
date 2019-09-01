#ifndef STAN_MATH_PRIM_SCAL_META_LENGTH_HPP
#define STAN_MATH_PRIM_SCAL_META_LENGTH_HPP

#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <cstdlib>

namespace stan {
/**
 * Returns the length of primitive scalar types
 * that are always of length 1.
 */
template <typename T, require_stan_scalar<T>...>
size_t length(T&& /*x*/) {
  return 1U;
}
}  // namespace stan
#endif
