#ifndef STAN_MATH_OPENCL_KERNELS_SCALAR_PARAMS_HPP
#define STAN_MATH_OPENCL_KERNELS_SCALAR_PARAMS_HPP
#ifdef STAN_OPENCL

namespace stan {
namespace math {
namespace opencl_kernels {

/** \ingroup opencl_kernels
 * Kernel source prefix for kernels that take their scalar parameters by
 * value. Kernels declare such parameters with `SCALAR_PARAM(name)` and read
 * them with `SCALAR_VALUE(name)`.
 */
static constexpr const char* scalar_params_by_value
    = "#define SCALAR_PARAM(name) const double name\n"
      "#define SCALAR_VALUE(name) name\n";

/** \ingroup opencl_kernels
 * Kernel source prefix for kernels that take their scalar parameters as
 * device scalars (`opencl::ScalarCl<double>`), which are buffers holding one
 * value.
 */
static constexpr const char* scalar_params_buffer
    = "#define SCALAR_PARAM(name) const __global double* name\n"
      "#define SCALAR_VALUE(name) name[0]\n";

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
