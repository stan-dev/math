#ifndef STAN_MATH_REV_FUN_EIGENVECTORS_HPP
#define STAN_MATH_REV_FUN_EIGENVECTORS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/fun/eigenvectors_sym.hpp>

namespace stan {
namespace math {

/**
 * Return the eigenvectors of the specified symmetric matrix.
 *
 * @tparam T type of input matrix.
 * @param m Specified matrix.
 * @return Eigenvectors of matrix.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto eigenvectors_sym(const T& m) {
  using return_t = return_var_matrix_t<T>;
  if (unlikely(m.size() == 0)) {
    return return_t(Eigen::MatrixXd(0, 0));
  }
  check_symmetric("eigenvectors_sym", "m", m);

  auto arena_m = to_arena(m);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(arena_m.val());
  arena_t<return_t> eigenvecs = solver.eigenvectors();
  auto eigenvals = to_arena(solver.eigenvalues());

  reverse_pass_callback([arena_m, eigenvals, eigenvecs]() mutable {
    const auto p = arena_m.val().cols();
    Eigen::MatrixXd f = (1
                         / (eigenvals.rowwise().replicate(p).transpose()
                            - eigenvals.rowwise().replicate(p))
                               .array());
    f.diagonal().setZero();
    arena_m.adj()
        += eigenvecs.val_op()
           * f.cwiseProduct(eigenvecs.val_op().transpose() * eigenvecs.adj_op())
           * eigenvecs.val_op().transpose();
  });

  return return_t(eigenvecs);
}

}  // namespace math
}  // namespace stan
#endif
