#ifndef STAN_MATH_REV_FUN_EIGENVECTORS_HPP
#define STAN_MATH_REV_FUN_EIGENVECTORS_HPP

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
 * Return the eigenvectors of the specified symmetric matrix.
 * <p>See <code>eigen_decompose()</code> for more information.
 * @tparam T type of input matrix.
 * @param m Specified matrix.
 * @return Eigenvectors of matrix.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto eigenvectors_sym(const T& m) {
  check_nonzero_size("eigenvalues_sym", "m", m);
  check_symmetric("eigenvalues_sym", "m", m);

  using return_t = return_var_matrix_t<T>;
  auto arena_m = to_arena(m);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(arena_m.val());
  arena_t<return_t> eigenvecs = solver.eigenvectors();
  auto eigenvals = to_arena(solver.eigenvalues());

  reverse_pass_callback([arena_m, eigenvals, eigenvecs]() mutable {
    const int p = arena_m.val().cols();
    auto eigenvecs_a = to_arena(eigenvecs.val());
    auto eigenvecs_a_adj = to_arena(eigenvecs.adj());

    matrix_d f(p, p);
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        f.coeffRef(j, i)
            = (i != j ? 1 / (eigenvals.coeff(i) - eigenvals.coeff(j)) : 0);

    arena_m.adj() += eigenvecs_a
                     * f.cwiseProduct(eigenvecs_a.transpose() * eigenvecs_a_adj)
                     * eigenvecs_a.transpose();
  });

  return return_t(eigenvecs);
}

}  // namespace math
}  // namespace stan
#endif
