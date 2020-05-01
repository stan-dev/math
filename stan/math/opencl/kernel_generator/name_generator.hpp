#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_NAME_GENERATOR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_NAME_GENERATOR_HPP
#ifdef STAN_OPENCL

#include <string>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */
/**
 * Unique name generator for variables used in generated kernels.
 */
class name_generator {
 public:
  /**
   * Generates a unique variable name.
   * @return variable name
   */
  inline std::string generate() { return "var" + std::to_string(++i_); }

 private:
  int i_{0};
};
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
