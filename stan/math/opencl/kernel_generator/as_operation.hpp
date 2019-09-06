#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Converts any valid kernel generator expression into an operation. Operations are returned unchanged.
 * @tparam T type of the input operation
 * @param a operation
 * @return operation
 */
template<typename T, typename = std::enable_if_t<std::is_base_of<operation_base, std::remove_reference_t<T>>::value>>
inline T&& as_operation(T&& a){
  return std::forward<T>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. Scalars are wrapped into \c scalar__.
 * @tparam T type of the input scalar
 * @param a scalar
 * @return \c scalar__ wrapping the input
 */
template<typename T, typename = enable_if_arithmetic<T>>
inline scalar__<T> as_operation(const T a){
  return scalar__<T>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. \c matrix_cl is wrapped into \c load__.
 * @tparam T \c matrix_cl
 * @param a \c matrix_cl
 * @return \c load__ wrapping the input
 */
template<typename T, typename = std::enable_if_t<std::is_base_of<matrix_cl<typename std::remove_reference_t<T>::type>,typename std::remove_reference_t<T>>::value>>
inline load__<T> as_operation(T&& a){
  return load__<T>(std::forward<T>(a));
}

/**
 * Type that results when converting any valid kernel generator expression into operation.
 */
template<typename T>
using as_operation_t = std::conditional_t<std::is_lvalue_reference<T>::value,
                                          decltype(as_operation(std::declval<T>())),
                                          std::remove_reference_t<decltype(as_operation(std::declval<T>()))>>;

}
}

#endif
#endif
