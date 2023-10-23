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
 * @tparam T type of input matrix.
 * @param m Specified matrix.
 * @return Eigenvalues of matrix.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto eigenvalues_sym(const T& m) {
  using return_t = return_var_matrix_t<Eigen::VectorXd, T>;
  if (unlikely(m.size() == 0)) {
    return return_t(Eigen::VectorXd(0));
  }
  check_symmetric("eigenvalues_sym", "m", m);

  auto arena_m = to_arena(m);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(arena_m.val());
  arena_t<return_t> eigenvals = solver.eigenvalues();
  auto eigenvecs = to_arena(solver.eigenvectors());

  reverse_pass_callback([eigenvals, arena_m, eigenvecs]() mutable {
    arena_m.adj()
        += eigenvecs * eigenvals.adj().asDiagonal() * eigenvecs.transpose();
  });

  return return_t(eigenvals);
}

}  // namespace math
}  // namespace stan
#endif
