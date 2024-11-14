#ifndef STAN_MATH_PRIM_META_IS_CONTAINER_HPP
#define STAN_MATH_PRIM_META_IS_CONTAINER_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_autodiff.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/base_type.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/**
 * Deduces whether type is eigen matrix or standard vector.
 * @tparam Container type to check
 */
template <typename Container>
using is_container = bool_constant<
    math::disjunction<is_eigen<Container>, is_std_vector<Container>>::value>;

/*! \ingroup general_types */
/*! \defgroup container_types container  */
/*! \addtogroup container_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_container */
/*! @tparam T the type to check */
template <typename T>
using require_container_t = require_t<is_container<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_container */
/*! @tparam T the type to check */
template <typename T>
using require_not_container_t = require_not_t<is_container<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_container */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_container_t
    = require_all_t<is_container<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_container */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_container_t
    = require_any_t<is_container<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_container */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_container_t
    = require_all_not_t<is_container<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup general_types */
/*! \defgroup container_types container  */
/*! \addtogroup container_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_container */
/*! and scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_container for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_container_st
    = require_t<container_type_check_base<is_container, scalar_type_t,
                                          TypeCheck, Check...>>;

/*! \brief Require type satisfies @ref is_container */
/*! and @ref base_type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_container for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_container_bt = require_t<
    container_type_check_base<is_container, base_type_t, TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_container */
/*! or scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_container for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_container_st
    = require_not_t<container_type_check_base<is_container, scalar_type_t,
                                              TypeCheck, Check...>>;

/*! \brief Require type satisfies @ref is_container */
/*! and holds a base type that satisfies @ref is_autodiff */
/*! @tparam T the type to check */
template <typename T>
using require_ad_container_t
    = require_all_t<stan::math::disjunction<is_eigen<T>, is_std_vector<T>>,
                    is_autodiff<base_type_t<T>>>;
/*! @} */

}  // namespace stan

#endif
