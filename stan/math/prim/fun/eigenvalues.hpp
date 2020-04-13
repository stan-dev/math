#ifndef STAN_MATH_PRIM_FUN_EIGENVALUES_HPP
#define STAN_MATH_PRIM_FUN_EIGENVALUES_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename T>
Eigen::Matrix<std::complex<T>, -1, 1> eigenvalues(
    const Eigen::Matrix<T, -1, -1>& m) {
  check_nonzero_size("eigenvalues", "m", m);
  check_square("eigenvalues", "m", m);

  Eigen::EigenSolver<Eigen::Matrix<T, -1, -1>> solver(m);
  return solver.eigenvalues();
}

}  // namespace math
}  // namespace stan
#endif
