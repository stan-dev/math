#ifndef STAN_MATH_PRIM_META_IS_STAN_SCALAR_HPP
#define STAN_MATH_PRIM_META_IS_STAN_SCALAR_HPP

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/**
 * Checks if decayed type is a var, fvar, or arithmetic
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_stan_scalar
    : bool_constant<math::disjunction<
          math::conjunction<is_var<std::decay_t<T>>,
                            std::is_arithmetic<value_type_t<T>>>,
          is_fvar<std::decay_t<T>>, std::is_arithmetic<std::decay_t<T>>,
          is_complex<std::decay_t<T>>>::value> {};

/*! \ingroup require_stan_scalar_real */
/*! \defgroup stan_scalar_types stan_scalar  */
/*! \addtogroup stan_scalar_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_stan_scalar */
/*! @tparam T the type to check */
template <typename T>
using require_stan_scalar_t = require_t<is_stan_scalar<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_stan_scalar */
/*! @tparam T the type to check */
template <typename T>
using require_not_stan_scalar_t
    = require_not_t<is_stan_scalar<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_stan_scalar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_stan_scalar_t
    = require_all_t<is_stan_scalar<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_stan_scalar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_stan_scalar_t
    = require_any_t<is_stan_scalar<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_stan_scalar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_stan_scalar_t
    = require_all_not_t<is_stan_scalar<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_stan_scalar
 */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_stan_scalar_t
    = require_any_not_t<is_stan_scalar<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \addtogroup stan_scalar_types */
/*! @{ */

/*! \brief Require value type does not satisfy @ref is_stan_scalar */
/*! @tparam T A type with a valid overload of @ref value_type available */
template <typename T>
using require_not_vt_stan_scalar
    = require_not_t<is_stan_scalar<value_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type satisfies @ref is_stan_scalar */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_st_stan_scalar
    = require_t<is_stan_scalar<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type does not satisfy @ref is_stan_scalar */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_not_st_stan_scalar
    = require_not_t<is_stan_scalar<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require all of the scalar types satisfy @ref is_stan_scalar */
/*! @tparam Types The types with a valid overload of @ref scalar_type available
 */
template <typename... Types>
using require_all_st_stan_scalar
    = require_all_t<is_stan_scalar<scalar_type_t<std::decay_t<Types>>>...>;
/*! @} */

}  // namespace stan

#endif
