#ifndef STAN_MATH_REV_FUN_EIGENVALUES_HPP
#define STAN_MATH_REV_FUN_EIGENVALUES_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
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
 *
 * @tparam T type of input matrix
 * @param m Input matrix.
 * @return Eigenvalues of matrix.
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto eigenvalues_sym(const T& m) {
  using ret_type = promote_scalar_t<var, Eigen::VectorXd>;

  check_square("eigenvalues_sym", "m", m);
  check_nonzero_size("eigenvalues_sym", "m", m);

  auto arena_m = to_arena(m);
  Eigen::MatrixXd m_val = value_of(arena_m);

  check_symmetric("eigenvalues_sym", "m", m_val);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m_val);
  arena_t<ret_type> eigenvalues = solver.eigenvalues();
  auto arena_eigenvectors_val = to_arena(solver.eigenvectors());

  reverse_pass_callback(
      [arena_m, arena_eigenvectors_val, eigenvalues]() mutable {
        arena_m.adj() += arena_eigenvectors_val * eigenvalues.adj().asDiagonal()
                         * arena_eigenvectors_val.transpose();
      });

  return ret_type(eigenvalues);
}

}  // namespace math
}  // namespace stan
#endif
