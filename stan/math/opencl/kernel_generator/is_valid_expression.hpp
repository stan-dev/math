#ifndef STAN_MATH_OPENCL_IS_USABLE_AS_OPERATION
#define STAN_MATH_OPENCL_IS_USABLE_AS_OPERATION
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 *Non-templated base of \c operation is needed for easy checking if something is
 *a subclass of \c operation.
 */
class operation_cl_base {};

/**
 * Determines whether a type is non-scalar type that is a valid kernel generator
 * expression.
 */
template <typename T, typename = void>
struct is_valid_expression_and_not_scalar
    : bool_constant<std::is_base_of<operation_cl_base,
                                    std::remove_reference_t<T>>::value> {};

template <typename T>
struct is_valid_expression_and_not_scalar<T, require_matrix_cl_t<T>>
    : std::true_type {};

/**
 * Determines whether a type is is a valid kernel generator expression. Valid
 * expressions are kernel generator operations, scalars and \c matric_cl and
 * references of these types.
 */
template <typename T>
struct is_valid_expression
    : bool_constant<is_valid_expression_and_not_scalar<T>::value
                    || std::is_arithmetic<std::remove_reference_t<T>>::value> {
};

/**
 * Enables a template if all given types are non-scalar types that are a
 * valid kernel generator expressions.
 */
template <typename... Types>
using require_all_valid_expressions_and_none_scalar_t
    = require_all_t<is_valid_expression_and_not_scalar<Types>...>;

/**
 * Enables a template if all given types are are a valid kernel generator
 * expressions.
 */
template <typename... Types>
using require_all_valid_expressions_t
    = require_all_t<is_valid_expression<Types>...>;

}  // namespace math
}  // namespace stan

#endif
#endif
