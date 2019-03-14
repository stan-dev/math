#ifndef STAN_MATH_OPENCL_CONSTANTS_HPP
#define STAN_MATH_OPENCL_CONSTANTS_HPP
#ifdef STAN_OPENCL
namespace stan {
namespace math {
enum class triangular_view_CL { LOWER = 0, UPPER = 1, ENTIRE = 2 };
enum class triangular_map_CL { UPPER_TO_LOWER = 0, LOWER_TO_UPPER = 1 };
}  // namespace math
}  // namespace stan
#endif
#endif
