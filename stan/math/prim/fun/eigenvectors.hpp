#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline auto eigenvectors(const EigMat& m) {
  using PlainMat = plain_type_t<EigMat>;
  const PlainMat& m_eval = m;
  check_nonzero_size("eigenvectors", "m", m_eval);
  check_square("eigenvectors", "m", m_eval);

  Eigen::EigenSolver<PlainMat> solver(m_eval);
  return solver.eigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
