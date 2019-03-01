#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_FVAR_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_FVAR_HPP

#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/or.hpp>

namespace stan {

/**
 * Metaprogram to calculate the base scalar return type resulting
 * from promoting all the scalar types of the template parameters.
 */
template <typename... T>
using contains_fvar = or_<is_fvar<typename scalar_type<T>::type>...>;

}  // namespace stan
#endif
