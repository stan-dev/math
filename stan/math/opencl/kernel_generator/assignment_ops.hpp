#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_ASSIGNMENT_OPS
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_ASSIGNMENT_OPS
#ifdef STAN_OPENCL
#include <stan/math/prim/meta/is_detected.hpp>

namespace stan {
namespace math {

/**
 * Ops that decide the type of assignment for LHS operations
 */
enum class assignment_ops_cl {equals, plus_equals, minus_equals, divide_equals};

/**
 * @param value A static constexpr const char* member for printing assignment ops
 */
template <assignment_ops_cl assign_op>
struct assignment_op_str;

template <>
struct assignment_op_str<assignment_ops_cl::equals> {
  static constexpr const char* value = " = ";
};

template <>
struct assignment_op_str<assignment_ops_cl::plus_equals> {
  static constexpr const char* value = " += ";
};

template <>
struct assignment_op_str<assignment_ops_cl::minus_equals> {
  static constexpr const char* value = " *= ";
};

template <>
struct assignment_op_str<assignment_ops_cl::divide_equals> {
  static constexpr const char* value = " /= ";
};


namespace internal {
template <typename, typename = void>
struct has_assignment_op_str : std::false_type {};

template <typename T>
struct has_assignment_op_str<T, void_t<decltype(T::assignment_op)>> : std::true_type {};

}  // namespace internal

/**
 * @tparam T A type that does not have an `assignment_op` static constexpr member type 
 * @return A constexpr const char* equal to `" = "`
 */
template <typename T, std::enable_if_t<!internal::has_assignment_op_str<std::decay_t<T>>::value>* = nullptr>
inline constexpr const char* assignment_op() noexcept {
  return " = ";
}

/**
 * @tparam T A type that has an `assignment_op` static constexpr member type
 * @return The types assignment op as a constexpr const char*
 */
template <typename T, std::enable_if_t<internal::has_assignment_op_str<T>::value>* = nullptr>
inline constexpr const char* assignment_op() noexcept {
  return assignment_op_str<std::decay_t<T>::assignment_op>::value;
}

}
}
#endif
#endif
