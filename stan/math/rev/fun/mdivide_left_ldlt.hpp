#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_alloc.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/core/chainable_object.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <memory>

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
  std::shared_ptr<Eigen::LDLT<Eigen::Matrix<double, R1, C1> > > ldltP_;
  Eigen::Matrix<double, R2, C2> C_;
};

/**
 * The vari for mdivide_left_ldlt(A, b) which handles the chain() call
 * for all elements of the result.  This vari follows the pattern
 * used in the other matrix operations where there is one "master"
 * vari whose value is never used and a large number of "slave" varis
 * whose chain() functions are never called because their adjoints are
 * set by the "master" vari.
 *
 * This class handles the var/var case.
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam R2 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C2 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 */
template <int R1, int C1, int R2, int C2>
class mdivide_left_ldlt_vv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefB_;
  vari **variRefC_;
  LDLT_factor<Eigen::Matrix<var, R1, C1>> A_;
  arena_t<Eigen::Matrix<double, R2, C2>> C_;

  mdivide_left_ldlt_vv_vari(const LDLT_factor<Eigen::Matrix<var, R1, C1>> &A,
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
        A_(A),
	C_(B.val()) {
    Eigen::Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    A_.solveInPlace(C_);
    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjB = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();

    A_.solveInPlace(adjB);

    A_.adj() -= adjB * C_.transpose();
    Eigen::Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

/**
 * The vari for mdivide_left_ldlt(A, b) which handles the chain() call
 * for all elements of the result.  This vari follows the pattern
 * used in the other matrix operations where there is one "master"
 * vari whose value is never used and a large number of "slave" varis
 * whose chain() functions are never called because their adjoints are
 * set by the "master" vari.
 *
 * This class handles the double/var case.
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam R2 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C2 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 */
template <int R1, int C1, int R2, int C2>
class mdivide_left_ldlt_dv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefB_;
  vari **variRefC_;
  chainable_object<LDLT_factor<Eigen::Matrix<double, R1, C1>>>* A_ptr_;
  arena_t<Eigen::Matrix<double, R2, C2>> C_;

  mdivide_left_ldlt_dv_vari(const LDLT_factor<Eigen::Matrix<double, R1, C1>> &A,
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
        A_ptr_(make_chainable_ptr(A)),
	C_(A_ptr_->get().solve(B.val())) {
    Eigen::Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjB = A_ptr_->get().solve(Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj());
    Eigen::Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

/**
 * The vari for mdivide_left_ldlt(A, b) which handles the chain() call
 * for all elements of the result.  This vari follows the pattern
 * used in the other matrix operations where there is one "master"
 * vari whose value is never used and a large number of "slave" varis
 * whose chain() functions are never called because their adjoints are
 * set by the "master" vari.
 *
 * This class handles the var/double case.
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam R2 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C2 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 */
template <int R1, int C1, int R2, int C2>
class mdivide_left_ldlt_vd_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefC_;
  LDLT_factor<Eigen::Matrix<var, R1, C1>> A_;
  arena_t<Eigen::Matrix<double, R2, C2>> C_;

  mdivide_left_ldlt_alloc<R1, C1, R2, C2> *alloc_;
  const LDLT_alloc<R1, C1> *alloc_ldlt_;

  mdivide_left_ldlt_vd_vari(const LDLT_factor<Eigen::Matrix<var, R1, C1>> &A,
                            const Eigen::Matrix<double, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
	A_(A),
	C_(B) {
    A_.solveInPlace(C_);
    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjC = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();

    A_.adj() -= A_.solve(adjC * C_.transpose());
  }
};
}  // namespace internal

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, typename EigMat,
          require_eigen_vt<is_var, EigMat> * = nullptr>
inline Eigen::Matrix<var, R1, EigMat::ColsAtCompileTime> mdivide_left_ldlt(
									   const LDLT_factor<Eigen::Matrix<var, R1, C1>> &A, const EigMat &b) {
  constexpr int R2 = EigMat::RowsAtCompileTime;
  constexpr int C2 = EigMat::ColsAtCompileTime;
  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);
  if (A.cols() == 0) {
    return {0, b.cols()};
  }

  auto *baseVari
      = new internal::mdivide_left_ldlt_vv_vari<R1, C1, R2, C2>(A, b);

  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, typename EigMat,
          require_eigen_vt<std::is_arithmetic, EigMat> * = nullptr>
inline Eigen::Matrix<var, R1, EigMat::ColsAtCompileTime> mdivide_left_ldlt(
									   const LDLT_factor<Eigen::Matrix<var, R1, C1>> &A, const EigMat &b) {
  constexpr int R2 = EigMat::RowsAtCompileTime;
  constexpr int C2 = EigMat::ColsAtCompileTime;
  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);
  if (A.cols() == 0) {
    return {0, b.cols()};
  }

  auto *baseVari
      = new internal::mdivide_left_ldlt_vd_vari<R1, C1, R2, C2>(A, b);

  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, typename EigMat,
          require_eigen_vt<is_var, EigMat> * = nullptr>
inline Eigen::Matrix<var, R1, EigMat::ColsAtCompileTime> mdivide_left_ldlt(
									   const LDLT_factor<Eigen::Matrix<double, R1, C1>> &A, const EigMat &b) {
  constexpr int R2 = EigMat::RowsAtCompileTime;
  constexpr int C2 = EigMat::ColsAtCompileTime;
  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);
  if (A.cols() == 0) {
    return {0, b.cols()};
  }

  auto *baseVari
      = new internal::mdivide_left_ldlt_dv_vari<R1, C1, R2, C2>(A, b);

  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

}  // namespace math
}  // namespace stan
#endif
