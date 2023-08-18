#ifndef STAN_MATH_REV_FUN_EIGENDECOMPOSE_HPP
#define STAN_MATH_REV_FUN_EIGENDECOMPOSE_HPP

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
 * Return the decomposition of the specified symmetric matrix.
 *
 * @tparam T type of input matrix.
 * @param m Specified matrix.
 * @return A tuple V,D where V is a matrix where the columns are the
 * eigenvectors of m, and D is a column vector of the eigenvalues of m.
 * The eigenvalues are in ascending order of magnitude, with the eigenvectors
 * provided in the same order.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto eigendecompose_sym(const T& m) {
  using eigval_return_t = return_var_matrix_t<Eigen::VectorXd, T>;
  using eigvec_return_t = return_var_matrix_t<T>;

  if (unlikely(m.size() == 0)) {
    return std::make_tuple(eigvec_return_t(Eigen::MatrixXd(0, 0)),
                           eigval_return_t(Eigen::VectorXd(0)));
  }
  check_symmetric("eigendecompose_sym", "m", m);

  auto arena_m = to_arena(m);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(arena_m.val());
  arena_t<eigval_return_t> eigenvals = solver.eigenvalues();
  arena_t<eigvec_return_t> eigenvecs = solver.eigenvectors();

  reverse_pass_callback([eigenvals, arena_m, eigenvecs]() mutable {
    // eigenvalue reverse calculation
    auto value_adj = eigenvecs.val_op() * eigenvals.adj().asDiagonal()
                     * eigenvecs.val_op().transpose();
    // eigenvector reverse calculation
    const auto p = arena_m.val().cols();
    Eigen::MatrixXd f
        = (1
           / (eigenvals.val_op().rowwise().replicate(p).transpose()
              - eigenvals.val_op().rowwise().replicate(p))
                 .array());
    f.diagonal().setZero();
    auto vector_adj
        = eigenvecs.val_op()
          * f.cwiseProduct(eigenvecs.val_op().transpose() * eigenvecs.adj_op())
          * eigenvecs.val_op().transpose();

    arena_m.adj() += value_adj + vector_adj;
  });

  return std::make_tuple(std::move(eigvec_return_t(eigenvecs)),
                         std::move(eigval_return_t(eigenvals)));
}

}  // namespace math
}  // namespace stan
#endif
