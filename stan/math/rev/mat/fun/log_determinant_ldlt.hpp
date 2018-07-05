#ifndef STAN_MATH_REV_MAT_FUN_LOG_DETERMINANT_LDLT_HPP
#define STAN_MATH_REV_MAT_FUN_LOG_DETERMINANT_LDLT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/LDLT_factor.hpp>

namespace stan {
namespace math {
namespace {

/**
 * Returns the log det of the matrix whose LDLT factorization is given
 * See The Matrix Cookbook's chapter on Derivatives of a Determinant
 * In this case, it is just the inverse of the underlying matrix
 * @param A, which is a LDLT_factor
 * @return ln(det(A))
 * @throws never
 */

template <int R, int C>
class log_det_ldlt_vari : public vari {
 public:
  explicit log_det_ldlt_vari(const LDLT_factor<var, R, C> &A)
      : vari(A.log_abs_det()),
        invA_mem_(ChainableStack::instance().memalloc_.alloc_array<double>(
            A.rows() * A.cols())),
        A_(A) {
    Eigen::Map<Eigen::Matrix<double, R, C>> invA(invA_mem_, A_.N_, A_.N_);

    invA.setIdentity();
    A_.solveInPlace(invA);
  }

  virtual void chain() {
    Eigen::Map<Eigen::Matrix<double, R, C>> invA(invA_mem_, A_.N_, A_.N_);

    for (size_t j = 0; j < A_.N_; j++) {
      for (size_t i = 0; i < A_.N_; i++) {
        A_.get_variA(i, j)->adj_ += adj_ * invA(i, j);
      }
    }
  }

  double *invA_mem_;
  LDLT_factor<var, R, C> A_;
};
}  // namespace

template <int R, int C>
var log_determinant_ldlt(LDLT_factor<var, R, C> &A) {
  return var(new log_det_ldlt_vari<R, C>(A));
}

}  // namespace math
}  // namespace stan
#endif
