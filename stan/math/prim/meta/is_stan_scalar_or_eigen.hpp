#ifndef STAN_MATH_PRIM_META_IS_STAN_SCALAR_OR_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_STAN_SCALAR_OR_EIGEN_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/** \ingroup type_trait
 * Extends std::true_type if all the provided types are either a Stan Scalar
 * type or a type inheriting from `EigenBase`.
 */
template <typename T>
using is_stan_scalar_or_eigen
    = bool_constant<is_stan_scalar<std::decay_t<T>>::value
                    || is_eigen<std::decay_t<T>>::value>;

//STAN_ADD_REQUIRE_UNARY(stan_scalar_or_eigen, is_stan_scalar_or_eigen, require_stan_scalar_real);
template <typename T>
using require_stan_scalar_or_eigen_t = require_t<is_stan_scalar_or_eigen<std::decay_t<T>>>;

template <typename T>
using require_not_stan_scalar_or_eigen_t
    = require_not_t<is_stan_scalar_or_eigen<std::decay_t<T>>>;

template <typename... Types>
using require_all_stan_scalar_or_eigen_t
    = require_all_t<is_stan_scalar_or_eigen<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_stan_scalar_or_eigen_t
    = require_any_t<is_stan_scalar_or_eigen<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_stan_scalar_or_eigen_t
    = require_all_not_t<is_stan_scalar_or_eigen<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_stan_scalar_or_eigen_t
    = require_any_not_t<is_stan_scalar_or_eigen<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(stan_scalar_or_eigen, is_stan_scalar_or_eigen, require_stan_scalar_real);
template <typename T>
using require_vt_stan_scalar_or_eigen
    = require_t<is_stan_scalar_or_eigen<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_stan_scalar_or_eigen
    = require_not_t<is_stan_scalar_or_eigen<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_stan_scalar_or_eigen
    = require_all_t<is_stan_scalar_or_eigen<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_stan_scalar_or_eigen
    = require_any_t<is_stan_scalar_or_eigen<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_stan_scalar_or_eigen
    = require_all_not_t<is_stan_scalar_or_eigen<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_stan_scalar_or_eigen
    = require_any_not_t<is_stan_scalar_or_eigen<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_stan_scalar_or_eigen
    = require_t<is_stan_scalar_or_eigen<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_stan_scalar_or_eigen
    = require_not_t<is_stan_scalar_or_eigen<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_stan_scalar_or_eigen
    = require_all_t<is_stan_scalar_or_eigen<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_stan_scalar_or_eigen
    = require_any_t<is_stan_scalar_or_eigen<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_stan_scalar_or_eigen
    = require_all_not_t<is_stan_scalar_or_eigen<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_stan_scalar_or_eigen
    = require_any_not_t<is_stan_scalar_or_eigen<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan
#endif
