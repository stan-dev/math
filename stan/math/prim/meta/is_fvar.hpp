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

//STAN_ADD_REQUIRE_UNARY(fvar, is_fvar, require_stan_scalar_real);
template <typename T>
using require_fvar_t = require_t<is_fvar<std::decay_t<T>>>;

template <typename T>
using require_not_fvar_t
    = require_not_t<is_fvar<std::decay_t<T>>>;

template <typename... Types>
using require_all_fvar_t
    = require_all_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_fvar_t
    = require_any_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_fvar_t
    = require_all_not_t<is_fvar<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_fvar_t
    = require_any_not_t<is_fvar<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(fvar, is_fvar, require_stan_scalar_real);
template <typename T>
using require_vt_fvar
    = require_t<is_fvar<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_fvar
    = require_not_t<is_fvar<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_fvar
    = require_all_t<is_fvar<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_fvar
    = require_any_t<is_fvar<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_fvar
    = require_all_not_t<is_fvar<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_fvar
    = require_any_not_t<is_fvar<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_fvar
    = require_t<is_fvar<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_fvar
    = require_not_t<is_fvar<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_fvar
    = require_all_t<is_fvar<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_fvar
    = require_any_t<is_fvar<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_fvar
    = require_all_not_t<is_fvar<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_fvar
    = require_any_not_t<is_fvar<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan
#endif
