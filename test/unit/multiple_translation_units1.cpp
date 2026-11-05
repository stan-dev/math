#include <test/unit/multiple_translation_units_includes.hpp>

inline stan::math::var function1() { return 0; }

#ifdef STAN_OPENCL
double opencl_cached_kernel1() {
  stan::math::matrix_cl<double> matrix = stan::math::constant(1.0, 2, 2);
  stan::math::matrix_cl<double> result = matrix * matrix;
  return stan::math::from_matrix_cl(result).sum();
}
#endif
