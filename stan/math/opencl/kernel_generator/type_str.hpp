#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_TYPE_STR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_TYPE_STR_HPP
#ifdef STAN_OPENCL

#include <string>

namespace stan {
namespace math {
/**
 * Determines a string name of a type. Unsupported types fail static assert.
 * @return name of the type
 */
template <typename T>
inline std::string type_str() {
  static_assert(sizeof(T) == -1, "Unsupported type in type_str");
  return "";
}

#define ADD_TYPE_TO_TYPE_STR(t)      \
  template <>                        \
  inline std::string type_str<t>() { \
    return #t;                       \
  }
ADD_TYPE_TO_TYPE_STR(double)
ADD_TYPE_TO_TYPE_STR(int)
ADD_TYPE_TO_TYPE_STR(char)
ADD_TYPE_TO_TYPE_STR(bool)
#undef ADD_TYPE_TO_TYPE_STR
}  // namespace math
}  // namespace stan
#endif
#endif
