#ifndef STAN_MATH_OPENCL_PRIM_EIGENVECTORS_SYM_HPP
#define STAN_MATH_OPENCL_PRIM_EIGENVECTORS_SYM_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/err/check_symmetric.hpp>
#include <stan/math/opencl/symmetric_eigensolver.hpp>

namespace stan {
namespace math {

inline matrix_cl<double> eigenvectors_sym(const matrix_cl<double>& m) {
  check_nonzero_size("eigenvectors_sym", "m", m);
  check_symmetric("eigenvalues_sym", "m", m);

  matrix_cl<double> eigenvalues, eigenvectors;
  symmetric_eigensolver(m, eigenvalues, eigenvectors);
  return eigenvectors;
}

}  // namespace math
}  // namespace stan
#endif
#endif
