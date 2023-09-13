#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_ASSIGNMENT_OPS
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_ASSIGNMENT_OPS
#ifdef STAN_OPENCL
#include <stan/math/prim/meta/is_detected.hpp>

namespace stan {
namespace math {

/**
 * Ops that decide the type of assignment for LHS operations
 */
enum class assign_op_cl {
  equals,
  plus_equals,
  minus_equals,
  divide_equals,
  multiply_equals
};

namespace internal {
/**
 * @param value A static constexpr const char* member for printing assignment
 * ops
 */
template <assign_op_cl assign_op>
struct assignment_op_str_impl;

template <>
struct assignment_op_str_impl<assign_op_cl::equals> {
  static constexpr const char* value = " = ";
};

template <>
struct assignment_op_str_impl<assign_op_cl::plus_equals> {
  static constexpr const char* value = " += ";
};

template <>
struct assignment_op_str_impl<assign_op_cl::minus_equals> {
  static constexpr const char* value = " -= ";
};

template <>
struct assignment_op_str_impl<assign_op_cl::divide_equals> {
  static constexpr const char* value = " /= ";
};

template <>
struct assignment_op_str_impl<assign_op_cl::multiply_equals> {
  static constexpr const char* value = " *= ";
};

template <typename, typename = void>
struct assignment_op_str : assignment_op_str_impl<assign_op_cl::equals> {};

template <typename T>
struct assignment_op_str<T, void_t<decltype(T::assignment_op)>>
    : assignment_op_str_impl<T::assignment_op> {};

}  // namespace internal

/**
 * @tparam T A type that has an `assignment_op` static constexpr member type
 * @return The types assignment op as a constexpr const char*
 */
template <typename T>
inline constexpr const char* assignment_op() noexcept {
  return internal::assignment_op_str<std::decay_t<T>>::value;
}

}  // namespace math
}  // namespace stan
#endif
#endif
