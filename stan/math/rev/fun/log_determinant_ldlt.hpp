#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_LDLT_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_alloc.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Returns the log det of the matrix whose LDLT factorization is given
 * See The Matrix Cookbook's chapter on Derivatives of a Determinant
 * In this case, it is just the inverse of the underlying matrix
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param A an LDLT_factor
 * @return ln(det(A))
 * @throws never
 */
template <int R, int C>
class log_det_ldlt_vari : public vari {
 public:
  explicit log_det_ldlt_vari(const LDLT_factor<var, R, C> &A)
      : vari(A.alloc_->log_abs_det()), alloc_ldlt_(A.alloc_) {}

  virtual void chain() {
    Eigen::Matrix<double, R, C> invA;

    // If we start computing Jacobians, this may be a bit inefficient
    invA.setIdentity(alloc_ldlt_->N_, alloc_ldlt_->N_);
    alloc_ldlt_->ldlt_.solveInPlace(invA);
    const_cast<matrix_vi &>(alloc_ldlt_->variA_).adj() += adj_ * invA;
  }
  const LDLT_alloc<R, C> *alloc_ldlt_;
};
}  // namespace internal

template <int R, int C>
var log_determinant_ldlt(LDLT_factor<var, R, C> &A) {
  if (A.rows() == 0) {
    return 0;
  }

  return var(new internal::log_det_ldlt_vari<R, C>(A));
}

}  // namespace math
}  // namespace stan
#endif
