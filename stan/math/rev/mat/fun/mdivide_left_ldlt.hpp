#ifndef STAN_MATH_REV_MAT_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_REV_MAT_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/LDLT_alloc.hpp>
#include <stan/math/rev/mat/fun/LDLT_factor.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>

namespace stan {
namespace math {
namespace internal {
template <int R1, int C1, int R2, int C2>
class mdivide_left_ldlt_alloc : public chainable_alloc {
 public:
  virtual ~mdivide_left_ldlt_alloc() {}

  /**
   * This share_ptr is used to prevent copying the LDLT factorizations
   * for mdivide_left_ldlt(ldltA, b) when ldltA is a LDLT_factor<double>.
   * The pointer is shared with the LDLT_factor<double> class.
   **/
  boost::shared_ptr<Eigen::LDLT<Eigen::Matrix<double, R1, C1> > > ldltP_;
  Eigen::Matrix<double, R2, C2> C_;
};

/**
 * The vari for mdivide_left_ldlt(A, b) which handles the chain() call
 * for all elements of the result.  This vari follows the pattern
 * used in the other matrix operations where there is one "master"
 * vari whose value is never used and a large number of "slave" varis
 * whose chain() functions are never called because their adjoints are
 * set by the "mater" vari.
 *
 * This class handles the var/var case.
 **/
template <int R1, int C1, int R2, int C2>
class mdivide_left_ldlt_vv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefB_;
  vari **variRefC_;
  mdivide_left_ldlt_alloc<R1, C1, R2, C2> *alloc_;
  const LDLT_alloc<R1, C1> *alloc_ldlt_;

  mdivide_left_ldlt_vv_vari(const LDLT_factor<var, R1, C1> &A,
                            const Eigen::Matrix<var, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        alloc_(new mdivide_left_ldlt_alloc<R1, C1, R2, C2>()),
        alloc_ldlt_(A.alloc_) {
    Eigen::Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    alloc_->C_ = B.val();
    alloc_ldlt_->ldlt_.solveInPlace(alloc_->C_);
    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = alloc_->C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjB = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();

    alloc_ldlt_->ldlt_.solveInPlace(adjB);

    const_cast<matrix_vi &>(alloc_ldlt_->variA_).adj()
        -= adjB * alloc_->C_.transpose();
    Eigen::Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

/**
 * The vari for mdivide_left_ldlt(A, b) which handles the chain() call
 * for all elements of the result.  This vari follows the pattern
 * used in the other matrix operations where there is one "master"
 * vari whose value is never used and a large number of "slave" varis
 * whose chain() functions are never called because their adjoints are
 * set by the "mater" vari.
 *
 * This class handles the double/var case.
 **/
template <int R1, int C1, int R2, int C2>
class mdivide_left_ldlt_dv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefB_;
  vari **variRefC_;
  mdivide_left_ldlt_alloc<R1, C1, R2, C2> *alloc_;

  mdivide_left_ldlt_dv_vari(const LDLT_factor<double, R1, C1> &A,
                            const Eigen::Matrix<var, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        alloc_(new mdivide_left_ldlt_alloc<R1, C1, R2, C2>()) {
    Eigen::Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    alloc_->C_ = B.val();
    alloc_->ldltP_ = A.ldltP_;
    alloc_->ldltP_->solveInPlace(alloc_->C_);
    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = alloc_->C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjB = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();
    alloc_->ldltP_->solveInPlace(adjB);
    Eigen::Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

/**
 * The vari for mdivide_left_ldlt(A, b) which handles the chain() call
 * for all elements of the result.  This vari follows the pattern
 * used in the other matrix operations where there is one "master"
 * vari whose value is never used and a large number of "slave" varis
 * whose chain() functions are never called because their adjoints are
 * set by the "mater" vari.
 *
 * This class handles the var/double case.
 **/
template <int R1, int C1, int R2, int C2>
class mdivide_left_ldlt_vd_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefC_;
  mdivide_left_ldlt_alloc<R1, C1, R2, C2> *alloc_;
  const LDLT_alloc<R1, C1> *alloc_ldlt_;

  mdivide_left_ldlt_vd_vari(const LDLT_factor<var, R1, C1> &A,
                            const Eigen::Matrix<double, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        alloc_(new mdivide_left_ldlt_alloc<R1, C1, R2, C2>()),
        alloc_ldlt_(A.alloc_) {
    alloc_->C_ = B;
    alloc_ldlt_->ldlt_.solveInPlace(alloc_->C_);
    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = alloc_->C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjC = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();

    const_cast<matrix_vi &>(alloc_ldlt_->variA_).adj()
        -= alloc_ldlt_->ldlt_.solve(adjC * alloc_->C_.transpose());
  }
};
}  // namespace internal

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<var, R1, C2> mdivide_left_ldlt(
    const LDLT_factor<var, R1, C1> &A, const Eigen::Matrix<var, R2, C2> &b) {
  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());

  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);

  internal::mdivide_left_ldlt_vv_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_ldlt_vv_vari<R1, C1, R2, C2>(A, b);

  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<var, R1, C2> mdivide_left_ldlt(
    const LDLT_factor<var, R1, C1> &A, const Eigen::Matrix<double, R2, C2> &b) {
  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());

  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);

  internal::mdivide_left_ldlt_vd_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_ldlt_vd_vari<R1, C1, R2, C2>(A, b);

  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<var, R1, C2> mdivide_left_ldlt(
    const LDLT_factor<double, R1, C1> &A, const Eigen::Matrix<var, R2, C2> &b) {
  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());

  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);

  internal::mdivide_left_ldlt_dv_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_ldlt_dv_vari<R1, C1, R2, C2>(A, b);

  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

}  // namespace math
}  // namespace stan
#endif
