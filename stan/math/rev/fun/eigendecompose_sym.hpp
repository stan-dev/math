#ifndef STAN_MATH_REV_FUN_EIGENDECOMPOSE_HPP
#define STAN_MATH_REV_FUN_EIGENDECOMPOSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/eigendecompose_sym.hpp>

#include <limits>

// W-40: relative eigenvalue gap below kappa * max(1, |w|_inf) * eps marks a
// numerically degenerate eigenvalue cluster for the reverse-mode adjoint.
#ifndef STAN_MATH_EIGEN_GAP_KAPPA
#define STAN_MATH_EIGEN_GAP_KAPPA 1e3
#endif

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
    auto value_adj = eigenvecs.val() * eigenvals.adj().asDiagonal()
                     * eigenvecs.val().transpose();
    // eigenvector reverse calculation
    const auto p = arena_m.val().cols();
    // W-40 cluster-aware adjoint: see rev/fun/eigenvectors_sym.hpp. Pairs
    // with an eigenvalue gap below STAN_MATH_EIGEN_GAP_KAPPA *
    // max(1, |w|_inf) * eps form a numerically degenerate cluster; their
    // 1/(w_j - w_i) coupling is replaced by zero (minimal-norm gauge).
    // Well-separated spectra take the original path VERBATIM below.
    const double eigengap_tau
        = STAN_MATH_EIGEN_GAP_KAPPA
          * std::max(1.0, eigenvals.val().cwiseAbs().maxCoeff())
          * std::numeric_limits<double>::epsilon();
    const bool has_degenerate_gap
        = p > 1
          && (eigenvals.val().tail(p - 1) - eigenvals.val().head(p - 1))
                     .minCoeff()
                 < eigengap_tau;
    Eigen::MatrixXd vector_adj;
    if (unlikely(has_degenerate_gap)) {
      Eigen::MatrixXd gaps
          = eigenvals.val().rowwise().replicate(p).transpose()
            - eigenvals.val().rowwise().replicate(p);
      Eigen::MatrixXd f = (gaps.array().abs() >= eigengap_tau)
                              .select(gaps.array().inverse(),
                                      Eigen::MatrixXd::Zero(p, p));
      vector_adj = eigenvecs.val()
                   * f.cwiseProduct(
                       eigenvecs.val().transpose() * eigenvecs.adj_op())
                   * eigenvecs.val().transpose();
    } else {
    Eigen::MatrixXd f = (1
                         / (eigenvals.val().rowwise().replicate(p).transpose()
                            - eigenvals.val().rowwise().replicate(p))
                               .array());
    f.diagonal().setZero();
    vector_adj = eigenvecs.val()
                 * f.cwiseProduct(
                     eigenvecs.val().transpose() * eigenvecs.adj_op())
                 * eigenvecs.val().transpose();
    }

    arena_m.adj() += value_adj + vector_adj;
  });

  return std::make_tuple(std::move(eigvec_return_t(eigenvecs)),
                         std::move(eigval_return_t(eigenvals)));
}

}  // namespace math
}  // namespace stan
#endif
