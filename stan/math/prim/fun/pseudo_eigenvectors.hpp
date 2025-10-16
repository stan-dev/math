#ifndef STAN_MATH_PRIM_FUN_PSEUDO_EIGENVECTORS_HPP
#define STAN_MATH_PRIM_FUN_PSEUDO_EIGENVECTORS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename T>
Eigen::Matrix<T, -1, -1> pseudo_eigenvectors(
    const Eigen::Matrix<T, -1, -1>& m) {
  check_nonzero_size("pseudo_eigenvectors", "m", m);
  check_square("pseudo_eigenvectors", "m", m);

  Eigen::EigenSolver<Eigen::Matrix<T, -1, -1>> solver(m);
  return solver.pseudoEigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
