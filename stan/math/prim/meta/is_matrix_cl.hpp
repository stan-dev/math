#ifndef STAN_MATH_PRIM_META_IS_MATRIX_CL_HPP
#define STAN_MATH_PRIM_META_IS_MATRIX_CL_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <typename T>
class arena_matrix_cl;

/**
 * Non-templated base class for `matrix_cl` simplifies checking if something is
 * matrix_cl.
 */
class matrix_cl_base {};

}  // namespace math

/** \ingroup type_traits
 * Checks if the decayed type of T is a matrix_cl.
 */
template <typename T>
struct is_matrix_cl
    : public std::is_base_of<math::matrix_cl_base, std::decay_t<T>> {};

/*! \ingroup matrix_cl_group */
/*! \defgroup matrix_cl_types matrix_cl  */
/*! \addtogroup matrix_cl_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_matrix_cl */
/*! @tparam T the type to check */
template <typename T>
using require_matrix_cl_t = require_t<is_matrix_cl<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_matrix_cl */
/*! @tparam T the type to check */
template <typename T>
using require_not_matrix_cl_t = require_not_t<is_matrix_cl<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_matrix_cl */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_matrix_cl_t
    = require_all_t<is_matrix_cl<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_matrix_cl */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_matrix_cl_t
    = require_any_t<is_matrix_cl<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_matrix_cl */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_matrix_cl_t
    = require_all_not_t<is_matrix_cl<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_matrix_cl */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_matrix_cl_t
    = require_any_not_t<is_matrix_cl<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup matrix_cl_group */
/*! \defgroup matrix_cl_types matrix_cl  */
/*! \addtogroup matrix_cl_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_matrix_cl */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_matrix_cl for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_matrix_cl_vt = require_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_matrix_cl or */
/*! value type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_matrix_cl for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_matrix_cl_vt = require_not_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_matrix_cl */
/*! and any of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_matrix_cl for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_matrix_cl_vt = require_any_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref is_matrix_cl
 */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_matrix_cl for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_matrix_cl_vt = require_any_not_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types satisfy @ref is_matrix_cl */
/*! and all of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_matrix_cl for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_matrix_cl_vt = require_all_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_matrix_cl */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_matrix_cl for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_matrix_cl_vt = require_all_not_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

/*! \brief Require type satisfies @ref is_matrix_cl */
/*! and scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_matrix_cl for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_matrix_cl_st
    = require_t<container_type_check_base<is_matrix_cl, scalar_type_t,
                                          TypeCheck, Check...>>;
/*! @} */

namespace internal {

template <typename T>
struct is_arena_matrix_cl_impl : public std::false_type {};

template <typename T>
struct is_arena_matrix_cl_impl<math::arena_matrix_cl<T>>
    : public std::true_type {};

}  // namespace internal

template <typename T>
struct is_arena_matrix_cl
    : public internal::is_arena_matrix_cl_impl<std::decay_t<T>> {};

}  // namespace stan
#endif
