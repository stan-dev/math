#ifndef STAN_MATH_PRIM_META_IS_AUTODIFF_HPP
#define STAN_MATH_PRIM_META_IS_AUTODIFF_HPP

#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <complex>
#include <type_traits>

namespace stan {

/**
 * Checks if decayed type is a var or fvar
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_autodiff
    : bool_constant<math::disjunction<is_var<std::decay_t<T>>,
                                      is_fvar<std::decay_t<T>>>::value> {};

//STAN_ADD_REQUIRE_UNARY(autodiff, is_autodiff, require_stan_scalar_real);
template <typename T>
using require_autodiff_t = require_t<is_autodiff<std::decay_t<T>>>;

template <typename T>
using require_not_autodiff_t
    = require_not_t<is_autodiff<std::decay_t<T>>>;

template <typename... Types>
using require_all_autodiff_t
    = require_all_t<is_autodiff<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_autodiff_t
    = require_any_t<is_autodiff<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_autodiff_t
    = require_all_not_t<is_autodiff<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_autodiff_t
    = require_any_not_t<is_autodiff<std::decay_t<Types>>...>;
  
//STAN_ADD_REQUIRE_UNARY_INNER(autodiff, is_autodiff, require_stan_scalar_real);
template <typename T>
using require_vt_autodiff
    = require_t<is_autodiff<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_autodiff
    = require_not_t<is_autodiff<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_autodiff
    = require_all_t<is_autodiff<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_autodiff
    = require_any_t<is_autodiff<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_autodiff
    = require_all_not_t<is_autodiff<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_autodiff
    = require_any_not_t<is_autodiff<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_autodiff
    = require_t<is_autodiff<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_autodiff
    = require_not_t<is_autodiff<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_autodiff
    = require_all_t<is_autodiff<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_autodiff
    = require_any_t<is_autodiff<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_autodiff
    = require_all_not_t<is_autodiff<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_autodiff
    = require_any_not_t<is_autodiff<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan

#endif
