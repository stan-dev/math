#ifndef STAN_MATH_REV_FUN_LDLT_ALLOC_HPP
#define STAN_MATH_REV_FUN_LDLT_ALLOC_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>

namespace stan {
namespace math {

/**
 * This object stores the actual (double typed) LDLT factorization of
 * an Eigen::Matrix<var> along with pointers to its vari's which allow the
 * *ldlt_ functions to save memory.  It is derived from a chainable_alloc
 * object so that it is allocated on the stack but does not have a chain()
 * function called.
 *
 * This class should only be instantiated as part of an LDLT_factor object
 * and is only used in *ldlt_ functions.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 **/
template <int R, int C>
class LDLT_alloc : public chainable_alloc {
 public:
  LDLT_alloc() : N_(0) {}
  explicit LDLT_alloc(const Eigen::Matrix<var, R, C> &A) : N_(0) { compute(A); }

  /**
   * Compute the LDLT factorization and store pointers to the
   * vari's of the matrix entries to be used when chain() is
   * called elsewhere.
   **/
  inline void compute(const Eigen::Matrix<var, R, C> &A) {
    N_ = A.rows();
    variA_ = A.vi();
    ldlt_.compute(A.val());
  }

  // Compute the log(abs(det(A))).  This is just a convenience function.
  inline double log_abs_det() const {
    return sum(log(ldlt_.vectorD().array()));
  }

  size_t N_;
  Eigen::LDLT<Eigen::Matrix<double, R, C> > ldlt_;
  Eigen::Matrix<vari *, R, C> variA_;
};

}  // namespace math
}  // namespace stan
#endif
