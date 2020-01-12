#ifndef STAN_MATH_PRIM_FUN_EIGENVALUES_SYM_HPP
#define STAN_MATH_PRIM_FUN_EIGENVALUES_SYM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the eigenvalues of the specified symmetric matrix
 * in descending order of magnitude.  This function is more
 * efficient than the general eigenvalues function for symmetric
 * matrices.
 * <p>See <code>eigen_decompose()</code> for more information.
 *
 * @tparam T type of elements in the matrix
 * @param m Specified matrix.
 * @return Eigenvalues of matrix.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> eigenvalues_sym(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  check_nonzero_size("eigenvalues_sym", "m", m);
  check_symmetric("eigenvalues_sym", "m", m);

  Eigen::SelfAdjointEigenSolver<
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
      solver(m, Eigen::EigenvaluesOnly);
  return solver.eigenvalues();
}

}  // namespace math
}  // namespace stan

#endif
