#ifndef STAN_MATH_PRIM_FUN_COMPLEX_SCHUR_DECOMPOSE_HPP
#define STAN_MATH_PRIM_FUN_COMPLEX_SCHUR_DECOMPOSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename T>
Eigen::Matrix<std::complex<T>, -1, -1> complex_schur_decompose(
    const Eigen::Matrix<T, -1, -1>& m) {
  check_nonzero_size("complex_schur_decompose", "m", m);
  check_square("complex_schur_decompose", "m", m);

  Eigen::ComplexSchur<Eigen::Matrix<T, -1, -1>> cs(m);
  return cs.matrixU();
}

}  // namespace math
}  // namespace stan
#endif
