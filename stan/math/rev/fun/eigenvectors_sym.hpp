#ifndef STAN_MATH_REV_FUN_EIGENVECTORS_HPP
#define STAN_MATH_REV_FUN_EIGENVECTORS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

namespace internal {
class eigenvectors_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols()
  double *A_;
  double *w_; // eigenvalues
  double *v_; // eigenvectors
  vari **vari_ref_A_;
  vari **vari_ref_w_;
  vari **vari_ref_v_;

  explicit eigenvectors_vari(const Eigen::Matrix<var, -1, -1> &A)
      : vari(0.0),
        M_(A.rows()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        w_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()))),
        v_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        vari_ref_A_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))),
        vari_ref_w_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()))),
        vari_ref_v_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))) {
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
    Map<matrix_vi>(vari_ref_v_, M_, M_)
        = vd.unaryExpr([](double x) { return new vari(x, false); });
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
std::cout << "HEHEHEHE" << std::endl;
    matrix_d adj_w = Map<matrix_vi>(vari_ref_w_, M_, 1).adj();
    matrix_d adj_v = Map<matrix_vi>(vari_ref_v_, M_, M_).adj();
    Map<matrix_d> w(w_, M_, M_);
    Map<matrix_d> v(v_, M_, M_);

    // f(i,j) = (i != j ? 1 / (e_i - e_j) : 0).
    matrix_d f(M_, M_);
    for (int i = 0; i < M_; i++)
      for (int j = 0; j < M_; j++)
        f.coeffRef(j ,i) = (i != j ? 1 / (w.coeff(i) - w.coeff(j)) : 0);
    std::cout << "F is\n" << f << std::endl;
    std::cout << "v.T * adj_v\n" << v.transpose() * adj_v << std::endl;
    std::cout << "F prod (v.T * adj_v)\n" << f.cwiseProduct(v.transpose() * adj_v) << std::endl;

    matrix_d diag_adj_w = adj_w.asDiagonal();

    std::cout << "diag_adj_w\n" << diag_adj_w << std::endl;

    std::cout << "adj_v\n" << adj_v << std::endl;

    matrix_d adjA =  v * f.cwiseProduct(v.transpose() * adj_v) * v.transpose();
    adjA += v * diag_adj_w * v.transpose();

    Map<matrix_vi>(vari_ref_A_, M_, M_).adj() += adjA;
  }
};
}  // namespace internal


/**
 * Return the eigenvectors of the specified symmetric matrix.
 * <p>See <code>eigen_decompose()</code> for more information.
 * @param m Specified matrix.
 * @return Eigenvectors of matrix.
 */
inline matrix_v eigenvectors_sym(const matrix_v & m) {
  check_symmetric("eigenvectors_sym", "m", m);
  matrix_v res(m.rows(), m.cols());
  internal::eigenvectors_vari *baseVari = new internal::eigenvectors_vari(m);
  res.vi() = Eigen::Map<matrix_vi>(baseVari->vari_ref_v_, res.rows(),
                                   res.cols());
  return res;
}

}  // namespace math
}  // namespace stan
#endif
