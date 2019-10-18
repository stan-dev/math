#ifndef STAN_MATH_REV_MAT_FUN_INVERSE_HPP
#define STAN_MATH_REV_MAT_FUN_INVERSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_nonempty.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

namespace stan {
namespace math {

namespace internal {
class inverse_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols()
  double *A_;
  double *A_inv_;
  vari **vari_ref_A_;
  vari **vari_ref_A_inv_;

  explicit inverse_vari(const Eigen::Matrix<var, -1, -1> &A)
      : vari(0.0),
        M_(A.rows()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        A_inv_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        vari_ref_A_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))),
        vari_ref_A_inv_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))) {
    using Eigen::Map;

    Map<matrix_d> Ad(A_, M_, M_);
    Map<matrix_d> A_inv_d(A_inv_, M_, M_);
    Ad = A.val();
    A_inv_d = Ad.inverse();

    Map<matrix_vi>(vari_ref_A_, M_, M_) = A.vi();
    Map<matrix_vi>(vari_ref_A_inv_, M_, M_)
        = A_inv_d.unaryExpr([](double x) { return new vari(x, false); });
  }

  /**
   * Reverse mode differentiation algorithm reference:
   *
   * Mike Giles. An extended collection of matrix derivative results for
   * forward and reverse mode AD.  Jan. 2008.
   *
   * Section 2.2.3 Inverse.
   */
  virtual void chain() {
    using Eigen::Map;

    matrix_d adj_A_inv = Map<matrix_vi>(vari_ref_A_inv_, M_, M_).adj();
    Map<matrix_d> A_inv_d(A_inv_, M_, M_);

    matrix_d adjA = A_inv_d.transpose() * adj_A_inv * A_inv_d.transpose();
    Map<matrix_vi>(vari_ref_A_, M_, M_).adj() -= adjA;
  }
};
}  // namespace internal

/**
 * Reverse mode specialization of calculating the inverse of the matrix.
 *
 * @param m Specified matrix.
 * @return Inverse of the matrix.
 */
inline matrix_v inverse(const matrix_v &m) {
  check_square("inverse", "m", m);
  check_nonempty("inverse", "m", m);
  matrix_v res(m.rows(), m.cols());
  internal::inverse_vari *baseVari = new internal::inverse_vari(m);
  res.vi() = Eigen::Map<matrix_vi>(baseVari->vari_ref_A_inv_, res.rows(),
                                   res.cols());
  return res;
}

}  // namespace math
}  // namespace stan
#endif
