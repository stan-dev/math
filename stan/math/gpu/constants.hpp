#ifndef STAN_MATH_GPU_CONSTANTS_HPP
#define STAN_MATH_GPU_CONSTANTS_HPP
#ifdef STAN_OPENCL
namespace stan {
namespace math {
namespace gpu {
const int Lower = 0;
const int Upper = 1;
const int Entire = 2;

const int LowerToUpper = 1;
const int UpperToLower = 0;
}  // namespace gpu
}  // namespace math
}  // namespace stan
#endif
#endif
