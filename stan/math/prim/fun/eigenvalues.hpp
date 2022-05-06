#ifndef STAN_MATH_PRIM_FUN_EIGENVALUES_HPP
#define STAN_MATH_PRIM_FUN_EIGENVALUES_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline auto eigenvalues(const EigMat& m) {
  check_nonzero_size("eigenvalues", "m", m);
  check_square("eigenvalues", "m", m);

  Eigen::EigenSolver<plain_type_t<EigMat>> solver(m);
  return solver.eigenvalues();
}

}  // namespace math
}  // namespace stan
#endif
