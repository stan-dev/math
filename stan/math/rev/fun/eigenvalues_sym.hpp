#ifndef STAN_MATH_REV_FUN_EIGENVALUES_HPP
#define STAN_MATH_REV_FUN_EIGENVALUES_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

namespace internal {
class eigenvalues_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols()
  double *A_;
  double *w_;  // eigenvalues
  double *v_;  // eigenvectors
  vari **vari_ref_A_;
  vari **vari_ref_w_;
  vari **vari_ref_v_;

  explicit eigenvalues_vari(const Eigen::Matrix<var, -1, -1> &A)
      : vari(0.0),
        M_(A.rows()),
        A_(ChainableStack::instance_->memalloc_.alloc_array<double>(
            A.rows() * A.cols())),
        w_(ChainableStack::instance_->memalloc_.alloc_array<double>(A.rows())),
        v_(ChainableStack::instance_->memalloc_.alloc_array<double>(
            A.rows() * A.cols())),
        vari_ref_A_(ChainableStack::instance_->memalloc_.alloc_array<vari *>(
            A.rows() * A.cols())),
        vari_ref_w_(ChainableStack::instance_->memalloc_.alloc_array<vari *>(
            A.rows())) {
    using Eigen::Map;

    Map<matrix_d> Ad(A_, M_, M_);
    Map<matrix_d> wd(w_, M_, 1);
    Map<matrix_d> vd(v_, M_, M_);
    Ad = A.val();
    Eigen::SelfAdjointEigenSolver<matrix_d> solver(Ad);
    wd = solver.eigenvalues();
    vd = solver.eigenvectors();

    Map<matrix_vi>(vari_ref_A_, M_, M_) = A.vi();
    Map<vector_vi>(vari_ref_w_, M_)
        = wd.unaryExpr([](double x) { return new vari(x, false); });
  }

  /**
   * Reverse mode differentiation algorithm reference:
   *
   * Mike Giles. An extended collection of matrix derivative results for
   * forward and reverse mode AD.  Jan. 2008.
   *
   * Section 3.1 Eigenvalues and eigenvectors.
   */
  virtual void chain() {
    using Eigen::Map;

    matrix_d adj_w = Map<matrix_vi>(vari_ref_w_, M_, 1).adj();
    Map<matrix_d> v(v_, M_, M_);

    matrix_d adjA = v * adj_w.asDiagonal() * v.transpose();

    Map<matrix_vi>(vari_ref_A_, M_, M_).adj() += adjA;
  }
};
}  // namespace internal

/**
 * Return the eigenvalues of the specified symmetric matrix.
 * <p>See <code>eigen_decompose()</code> for more information.
 * @param m Specified matrix.
 * @return Eigenvalues of matrix.
 */
inline vector_v eigenvalues_sym(const matrix_v &m) {
  matrix_d m_eval(value_of_rec(m));
  check_nonzero_size("eigenvalues_sym", "m", m_eval);
  check_symmetric("eigenvalues_sym", "m", m_eval);
  vector_v res(m.rows());
  internal::eigenvalues_vari *baseVari = new internal::eigenvalues_vari(m);
  res.vi()
      = Eigen::Map<vector_vi>(baseVari->vari_ref_w_, res.rows(), res.cols());
  return res;
}

}  // namespace math
}  // namespace stan
#endif
