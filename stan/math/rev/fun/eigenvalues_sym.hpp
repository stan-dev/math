#ifndef STAN_MATH_REV_FUN_EIGENVALUES_HPP
#define STAN_MATH_REV_FUN_EIGENVALUES_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Return the eigenvalues of the specified symmetric matrix.
 * <p>See <code>eigen_decompose()</code> for more information.
 * @param m Specified matrix.
 * @return Eigenvalues of matrix.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto eigenvalues_sym(const T& m) {
  check_nonzero_size("eigenvalues_sym", "m", m);
  check_symmetric("eigenvalues_sym", "m", m);

  using return_t = return_var_matrix_t<T>;
  auto arena_m = to_arena(m);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(arena_m.val());
  arena_t<return_t> eigenvals = solver.eigenvalues();
  auto eigenvecs = to_arena(solver.eigenvectors());
  auto eigenvals_adj = eigenvals.adj();

  reverse_pass_callback([eigenvals, arena_m, eigenvecs, eigenvals_adj]() mutable {
    arena_m.adj() +=
        eigenvecs * eigenvals_adj.asDiagonal() * eigenvecs.transpose();
  });

  return return_t(eigenvals);
}

}  // namespace math
}  // namespace stan
#endif


