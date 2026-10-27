#include <test/unit/multiple_translation_units_includes.hpp>

inline stan::math::var function2() { return 0; }

#ifdef STAN_OPENCL
double opencl_cached_kernel2() {
  stan::math::matrix_cl<double> matrix = stan::math::constant(2.0, 2, 2);
  stan::math::matrix_cl<double> result = matrix + matrix;
  return stan::math::from_matrix_cl(result).sum();
}

void select_opencl_device2() {
  stan::math::opencl_context.select_device(OPENCL_PLATFORM_ID,
                                           OPENCL_DEVICE_ID);
}
#endif
