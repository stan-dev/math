#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename EigMat>
Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>
eigenvectors_sym(const EigMat& m) {
  const plain_type_t<EigMat>& m_eval = m;
  check_nonzero_size("eigenvectors_sym", "m", m_eval);
  check_symmetric("eigenvalues_sym", "m", m_eval);

  Eigen::SelfAdjointEigenSolver<plain_type_t<EigMat>> solver(m_eval);
  return solver.eigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
