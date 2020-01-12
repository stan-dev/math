#ifndef STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP
#define STAN_MATH_PRIM_FUN_EIGENVECTORS_SYM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenvectors_sym(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  check_nonzero_size("eigenvectors_sym", "m", m);
  check_symmetric("eigenvalues_sym", "m", m);

  Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
      solver(m);
  return solver.eigenvectors();
}

}  // namespace math
}  // namespace stan
#endif
