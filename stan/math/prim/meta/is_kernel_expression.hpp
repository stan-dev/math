#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_IS_KERNEL_EXPRESSION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_IS_KERNEL_EXPRESSION_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_matrix_cl.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Non-templated base of `operation_cl` is needed for easy checking if something
 * is a subclass of `operation_cl`.
 */
class operation_cl_base {};
/**
 * Non-templated base of `operation_cl_lhs` is needed for easy checking if
 * something is a subclass of `operation_cl_lhs`.
 */
class operation_cl_lhs_base {};

/**
 * Determines whether a type is non-scalar type that is a valid kernel generator
 * expression.
 */
template <typename T, typename = void>
struct is_kernel_expression_and_not_scalar
    : bool_constant<std::is_base_of<operation_cl_base,
                                    std::remove_reference_t<T>>::value> {};
template <typename T>
struct is_kernel_expression_and_not_scalar<T, require_matrix_cl_t<T>>
    : std::true_type {};

/**
 * Determines whether a type is is a valid kernel generator expression. Valid
 * expressions are kernel generator operations, scalars and \c matrix_cl and
 * references of these types.
 */
template <typename T>
struct is_kernel_expression
    : bool_constant<is_kernel_expression_and_not_scalar<T>::value
                    || std::is_arithmetic<std::remove_reference_t<T>>::value> {
};

/**
 * Enables a template if all given types are non-scalar types that are a
 * valid kernel generator expressions.
 */
template <typename... Types>
using require_all_kernel_expressions_and_none_scalar_t
    = require_all_t<is_kernel_expression_and_not_scalar<Types>...>;

/**
 * Enables a template if all given types are are a valid kernel generator
 * expressions.
 */
template <typename... Types>
using require_all_kernel_expressions_t
    = require_all_t<is_kernel_expression<Types>...>;

/**
 * Determines whether a type is an assignable kernel generator
 * expression.
 */
template <typename T, typename = void>
struct is_kernel_expression_lhs
    : bool_constant<std::is_base_of<operation_cl_lhs_base,
                                    std::remove_reference_t<T>>::value> {};
template <typename T>
struct is_kernel_expression_lhs<T, require_matrix_cl_t<T>> : std::true_type {};

/**
 * Determines whether a type is a var containing a kernel generator expression.
 */
template <typename T>
struct is_rev_kernel_expression
    : math::conjunction<is_var<T>, is_kernel_expression<value_type_t<T>>> {};

/**
 * Determines whether a type is either a kernel generator
 * expression or a var containing a kernel generator expression.
 */
template <typename T>
struct is_prim_or_rev_kernel_expression
    : math::disjunction<is_kernel_expression<T>, is_rev_kernel_expression<T>> {
};

/**
 * Determines whether a type is either a non-scalar kernel generator
 * expression or a var containing a non-scalar kernel generator expression.
 */
template <typename T>
struct is_nonscalar_prim_or_rev_kernel_expression
    : math::disjunction<
          is_kernel_expression_and_not_scalar<T>,
          math::conjunction<is_var<T>, is_kernel_expression_and_not_scalar<
                                           value_type_t<T>>>> {};

/** @}*/

/*! \ingroup opencl_kernel_generator */
/*! \defgroup kernel_expression_lhs_types kernel_expression_lhs  */
/*! \addtogroup kernel_expression_lhs_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_kernel_expression_lhs */
/*! @tparam T the type to check */
template <typename T>
using require_kernel_expression_lhs_t
    = require_t<is_kernel_expression_lhs<std::decay_t<T>>>;
/*! @} */

/*! \ingroup opencl_kernel_generator */
/*! \defgroup rev_kernel_expression_types rev_kernel_expression  */
/*! \addtogroup rev_kernel_expression_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_rev_kernel_expression */
/*! @tparam T the type to check */
template <typename T>
using require_rev_kernel_expression_t
    = require_t<is_rev_kernel_expression<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_rev_kernel_expression */
/*! @tparam T the type to check */
template <typename T>
using require_not_rev_kernel_expression_t
    = require_not_t<is_rev_kernel_expression<std::decay_t<T>>>;
/*! @} */

/*! \ingroup opencl_kernel_generator */
/*! \defgroup prim_or_rev_kernel_expression_types prim_or_rev_kernel_expression
 */
/*! \addtogroup prim_or_rev_kernel_expression_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_prim_or_rev_kernel_expression */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_prim_or_rev_kernel_expression_t
    = require_all_t<is_prim_or_rev_kernel_expression<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup opencl_kernel_generator */
/*! \defgroup nonscalar_prim_or_rev_kernel_expression_types
 * nonscalar_prim_or_rev_kernel_expression  */
/*! \addtogroup nonscalar_prim_or_rev_kernel_expression_types */
/*! @{ */

/*! \brief Require type satisfies @ref
 * is_nonscalar_prim_or_rev_kernel_expression */
/*! @tparam Types The types that are checked */
template <typename T>
using require_nonscalar_prim_or_rev_kernel_expression_t
    = require_t<is_nonscalar_prim_or_rev_kernel_expression<std::decay_t<T>>>;

/*! \brief Require type does not satisfy
 * @ref is_nonscalar_prim_or_rev_kernel_expression */
/*! @tparam T the type to check */
template <typename T>
using require_not_nonscalar_prim_or_rev_kernel_expression_t = require_not_t<
    is_nonscalar_prim_or_rev_kernel_expression<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy
 * @ref is_nonscalar_prim_or_rev_kernel_expression */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_nonscalar_prim_or_rev_kernel_expression_t = require_all_t<
    is_nonscalar_prim_or_rev_kernel_expression<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy
 * @ref is_nonscalar_prim_or_rev_kernel_expression */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_nonscalar_prim_or_rev_kernel_expression_t = require_any_t<
    is_nonscalar_prim_or_rev_kernel_expression<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy
 * @ref is_nonscalar_prim_or_rev_kernel_expression */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_nonscalar_prim_or_rev_kernel_expression_t
    = require_all_not_t<
        is_nonscalar_prim_or_rev_kernel_expression<std::decay_t<Types>>...>;
/*! @} */
}  // namespace stan

#endif
