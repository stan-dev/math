#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_FVAR_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_FVAR_HPP

#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/or.hpp>

namespace stan {

/**
 * Metaprogram to determine if any of the
 * provided types is a fvar.
 */
template <typename... T>
using contains_fvar = math::or_<is_fvar<typename scalar_type<T>::type>...>;

}  // namespace stan
#endif
