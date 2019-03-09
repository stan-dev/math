#ifndef STAN_MATH_PRIM_SCAL_META_CONTAINS_FVAR_HPP
#define STAN_MATH_PRIM_SCAL_META_CONTAINS_FVAR_HPP

#include <stan/math/prim/scal/meta/is_fvar.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

namespace stan {

/**
 * Defines a public enum named value which is defined to be true (1)
 * if any of the template parameters includes a fvar as their base scalar and
 * false (0) otherwise.
 */
template <typename... T>
using contains_fvar
    = math::disjunction<is_fvar<typename scalar_type<T>::type>...>;

}  // namespace stan
#endif
