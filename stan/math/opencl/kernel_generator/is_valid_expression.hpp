#ifndef STAN_MATH_OPENCL_IS_USABLE_AS_OPERATION
#define STAN_MATH_OPENCL_IS_USABLE_AS_OPERATION
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Determines whether a type is non-scalar type that is a valid kernel generator expression.
 */
template<typename T>
struct is_valid_expression_and_not_scalar{
  enum{ value = std::is_base_of<operation_base,T>::value};
};

template<typename T>
struct is_valid_expression_and_not_scalar<matrix_cl < T>> : std::true_type{};

template<typename T>
struct is_valid_expression_and_not_scalar<T&>{
    enum{ value = is_valid_expression_and_not_scalar<T>::value};
};

template<typename T>
struct is_valid_expression_and_not_scalar<T&&>{
    enum{ value = is_valid_expression_and_not_scalar<T>::value};
};

/**
 * Determines whether a type is is a valid kernel generator expression. Valid expressions are kernel generator operations, scalars and \c matric_cl and references of these types.
 */
template<typename T>
struct is_valid_expression{
    enum{ value = is_valid_expression_and_not_scalar<T>::value || std::is_arithmetic<std::remove_reference_t<T>>::value};
};

/**
 * Disables a template if not all given types are non-scalar types that are a valid kernel generator expressions.
 */
template<typename... Types>
using enable_if_all_valid_expressions_and_none_scalar = typename std::enable_if_t<math::conjunction<is_valid_expression_and_not_scalar < Types>...>::value>;

/**
 * Disables a template if not all given types are are a valid kernel generator expressions.
 */
template<typename... Types>
using enable_if_all_valid_expressions = typename std::enable_if_t<math::conjunction<is_valid_expression<Types>...>::value>;

}
}


#endif
#endif
