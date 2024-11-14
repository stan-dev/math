#ifndef STAN_MATH_PRIM_FUN_PSEUDO_EIGENVECTORS_HPP
#define STAN_MATH_PRIM_FUN_PSEUDO_EIGENVECTORS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_matrix_base_t<EigMat>* = nullptr>
inline Eigen::Matrix<scalar_type_t<EigMat>, -1, -1> pseudo_eigenvectors(
    EigMat&& m) {
  check_nonzero_size("pseudo_eigenvectors", "m", m);
  check_square("pseudo_eigenvectors", "m", m);
  Eigen::EigenSolver<Eigen::Matrix<scalar_type_t<EigMat>, -1, -1>> solver(
      to_ref(std::forward<EigMat>(m)));
  return solver.pseudoEigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
