#ifndef STAN_MATH_PRIM_META_IS_REAL_HPP
#define STAN_MATH_PRIM_META_IS_REAL_HPP

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
 * Checks if decayed type is a var, fvar, or arithmetic. This
 * is equivalent to @ref is_stan_scalar, except for excluding std::complex
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_real : bool_constant<math::disjunction<
                     math::conjunction<is_var<std::decay_t<T>>,
                                       std::is_arithmetic<value_type_t<T>>>,
                     is_fvar<std::decay_t<T>>,
                     std::is_arithmetic<std::decay_t<T>>>::value> {};

template <typename T>
constexpr bool is_real_v = is_real<T>::value;

/*! \ingroup require_real_real */
/*! \defgroup real_types real  */
/*! \addtogroup real_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_real */
/*! @tparam T the type to check */
template <typename T>
using require_real_t = require_t<is_real<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_real */
/*! @tparam T the type to check */
template <typename T>
using require_not_real_t = require_not_t<is_real<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_real */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_real_t = require_all_t<is_real<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_real */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_real_t = require_any_t<is_real<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_real */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_real_t
    = require_all_not_t<is_real<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_real
 */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_real_t
    = require_any_not_t<is_real<std::decay_t<Types>>...>;
/*! @} */

}  // namespace stan

#endif
