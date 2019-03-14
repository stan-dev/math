#ifndef STAN_MATH_OPENCL_CONSTANTS_HPP
#define STAN_MATH_OPENCL_CONSTANTS_HPP
#ifdef STAN_OPENCL
namespace stan {
namespace math {
enum class triangular_view_CL { LOWER = 0, UPPER = 1, ENTIRE = 2 };
enum class triangular_map_CL { UPPER_TO_LOWER = 0, LOWER_TO_UPPER = 1 };

namespace opencl_kernels {
/**
 * Returns the enum string for OpenCL kernels
 *
 * @return enum string
 */
inline std::string get_opencl_enum_string() {
  std::ostringstream enums;
  enums << "typedef enum {\n";
  enums << "    LOWER=" << static_cast<int>(triangular_view_CL::LOWER) << ",\n";
  enums << "    UPPER=" << static_cast<int>(triangular_view_CL::UPPER) << ",\n";
  enums << "    ENTIRE=" << static_cast<int>(triangular_view_CL::ENTIRE)
        << "\n";
  enums << "} triangular_view;";
  enums << "typedef enum {\n";
  enums << "    UPPER_TO_LOWER="
        << static_cast<int>(triangular_map_CL::UPPER_TO_LOWER) << ",\n";
  enums << "    LOWER_TO_UPPER="
        << static_cast<int>(triangular_map_CL::LOWER_TO_UPPER) << "\n";
  enums << "} triangular_map;";
  return enums.str();
}
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
