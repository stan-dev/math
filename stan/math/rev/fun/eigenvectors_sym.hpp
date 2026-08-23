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

#include <limits>

// W-40: relative eigenvalue gap below kappa * max(1, |w|_inf) * eps marks a
// numerically degenerate eigenvalue cluster for the reverse-mode adjoint.
#ifndef STAN_MATH_EIGEN_GAP_KAPPA
#define STAN_MATH_EIGEN_GAP_KAPPA 1e3
#endif

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
    // Eigenvalue pairs with a gap below STAN_MATH_EIGEN_GAP_KAPPA *
    // max(1, |w|_inf) * eps are treated as one numerically degenerate
    // cluster: the individual eigenvectors are not identifiable at machine
    // resolution (any orthogonal basis of the cluster subspace is a valid
    // answer), and the standard adjoint's 1/(w_j - w_i) coupling amplifies
    // rounding-level differences to O(1/eps). For such pairs the coupling
    // is replaced by its minimal-norm (zero) value -- the cluster-gauge
    // choice; cross-cluster pairs keep the standard formula. Well-separated
    // spectra (min adjacent gap >= threshold) take the original code path
    // below VERBATIM, so results there are bit-identical to the unguarded
    // implementation.
    const double eigengap_tau = STAN_MATH_EIGEN_GAP_KAPPA
                                * std::max(1.0, eigenvals.cwiseAbs().maxCoeff())
                                * std::numeric_limits<double>::epsilon();
    const bool has_degenerate_gap
        = p > 1 && (eigenvals.tail(p - 1) - eigenvals.head(p - 1)).minCoeff()
                       < eigengap_tau;
    if (unlikely(has_degenerate_gap)) {
      Eigen::MatrixXd gaps = eigenvals.rowwise().replicate(p).transpose()
                             - eigenvals.rowwise().replicate(p);
      Eigen::MatrixXd f = (gaps.array().abs() >= eigengap_tau)
                              .select(gaps.array().inverse(),
                                      Eigen::MatrixXd::Zero(p, p));
      arena_m.adj()
          += eigenvecs.val()
             * f.cwiseProduct(eigenvecs.val().transpose() * eigenvecs.adj_op())
             * eigenvecs.val().transpose();
    } else {
    Eigen::MatrixXd f = (1
                         / (eigenvals.rowwise().replicate(p).transpose()
                            - eigenvals.rowwise().replicate(p))
                               .array());
    f.diagonal().setZero();
    arena_m.adj()
        += eigenvecs.val()
           * f.cwiseProduct(eigenvecs.val().transpose() * eigenvecs.adj_op())
           * eigenvecs.val().transpose();
    }
  });

  return return_t(eigenvecs);
}

}  // namespace math
}  // namespace stan
#endif
