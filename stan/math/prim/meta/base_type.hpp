#ifndef STAN_MATH_PRIM_META_BASE_TYPE_HPP
#define STAN_MATH_PRIM_META_BASE_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/**
 * Metaprogram structure to determine the base base type of a template
 * argument.  Qualifiers `const` and `volatile` are removed from all
 * types as are references.
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T type of non-container
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct base_type {
  using type = std::decay_t<T>;
};

template <typename T>
using base_type_t = typename base_type<T>::type;

/**
 * Specialization of base_type for vector to recursively return the inner
 * base type.
 *
 * @tparam T type of standard vector
 * @ingroup type_trait
 */
template <typename T>
struct base_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = base_type_t<typename std::decay_t<T>::value_type>;
};

/**
 * Template metaprogram defining the base base type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type of matrix
 * @ingroup type_trait
 */
template <typename T>
struct base_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = base_type_t<typename std::decay_t<T>::Scalar>;
};

/*! \ingroup require_stan_scalar_real */
/*! \addtogroup var_types */
/*! @{ */

/*! \brief Require any of the base types satisfy @ref is_var */
/*! Unlike @ref require_any_st_var this also matches `std::complex<var>`. */
/*! Defined here because is_var.hpp cannot include base_type.hpp. */
/*! @tparam Types The types with a valid overload of @ref base_type available
 */
template <typename... Types>
using require_any_bt_var
    = require_any_t<is_var<base_type_t<std::decay_t<Types>>>...>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \addtogroup arithmetic_types */
/*! @{ */

/*! \brief Require all of the base types satisfy std::is_arithmetic */
/*! Unlike @ref require_all_st_arithmetic this also matches */
/*! `std::complex<double>`. */
/*! @tparam Types The types with a valid overload of @ref base_type available
 */
template <typename... Types>
using require_all_bt_arithmetic
    = require_all_t<std::is_arithmetic<base_type_t<std::decay_t<Types>>>...>;
/*! @} */

}  // namespace stan
#endif
