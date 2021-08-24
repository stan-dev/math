#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename T>
inline auto eigenvectors(const Eigen::Matrix<T, -1, -1>& m) {
  check_nonzero_size("eigenvectors", "m", m);
  check_square("eigenvectors", "m", m);

  Eigen::EigenSolver<Eigen::Matrix<T, -1, -1>> solver(m);
  return solver.eigenvectors().eval();
}

}  // namespace math
}  // namespace stan
#endif
