#ifndef STAN_MATH_REV_MAT_FUN_MDIVIDE_LEFT_SPD_HPP
#define STAN_MATH_REV_MAT_FUN_MDIVIDE_LEFT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_alloc : public chainable_alloc {
 public:
  virtual ~mdivide_left_spd_alloc() {}

  Eigen::LLT<Eigen::Matrix<double, R1, C1> > llt_;
  Eigen::Matrix<double, R2, C2> C_;
};

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_vv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefA_;
  vari **variRefB_;
  vari **variRefC_;
  mdivide_left_spd_alloc<R1, C1, R2, C2> *alloc_;

  mdivide_left_spd_vv_vari(const Eigen::Matrix<var, R1, C1> &A,
                           const Eigen::Matrix<var, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        variRefA_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        alloc_(new mdivide_left_spd_alloc<R1, C1, R2, C2>()) {
    Eigen::Map<matrix_vi>(variRefA_, M_, M_) = A.vi();
    Eigen::Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    alloc_->C_ = B.val();
    alloc_->llt_ = A.val().llt();
    alloc_->llt_.solveInPlace(alloc_->C_);

    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = alloc_->C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjB = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();
    alloc_->llt_.solveInPlace(adjB);
    Eigen::Map<matrix_vi>(variRefA_, M_, M_).adj()
        -= adjB * alloc_->C_.transpose();
    Eigen::Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_dv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefB_;
  vari **variRefC_;
  mdivide_left_spd_alloc<R1, C1, R2, C2> *alloc_;

  mdivide_left_spd_dv_vari(const Eigen::Matrix<double, R1, C1> &A,
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
        alloc_(new mdivide_left_spd_alloc<R1, C1, R2, C2>()) {
    alloc_->C_ = B.val();
    Eigen::Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    alloc_->llt_ = A.llt();
    alloc_->llt_.solveInPlace(alloc_->C_);

    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = alloc_->C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjB = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();
    alloc_->llt_.solveInPlace(adjB);
    Eigen::Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_vd_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  vari **variRefA_;
  vari **variRefC_;
  mdivide_left_spd_alloc<R1, C1, R2, C2> *alloc_;

  mdivide_left_spd_vd_vari(const Eigen::Matrix<var, R1, C1> &A,
                           const Eigen::Matrix<double, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        variRefA_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        alloc_(new mdivide_left_spd_alloc<R1, C1, R2, C2>()) {
    Eigen::Map<matrix_vi>(variRefA_, M_, M_) = A.vi();
    alloc_->llt_ = A.val().llt();
    alloc_->C_ = alloc_->llt_.solve(B);

    Eigen::Map<matrix_vi>(variRefC_, M_, N_)
        = alloc_->C_.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    matrix_d adjC = Eigen::Map<matrix_vi>(variRefC_, M_, N_).adj();
    Eigen::Map<matrix_vi>(variRefA_, M_, M_).adj()
        -= alloc_->llt_.solve(adjC * alloc_->C_.transpose());
  }
};
}  // namespace internal

template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<var, R1, C2> mdivide_left_spd(
    const Eigen::Matrix<var, R1, C1> &A, const Eigen::Matrix<var, R2, C2> &b) {
  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());

  check_square("mdivide_left_spd", "A", A);
  check_multiplicable("mdivide_left_spd", "A", A, "b", b);

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  internal::mdivide_left_spd_vv_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_spd_vv_vari<R1, C1, R2, C2>(A, b);

  res.vi() = Eigen::Map<matrix_vi>(&baseVari->variRefC_[0], b.rows(), b.cols());
  return res;
}

template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<var, R1, C2> mdivide_left_spd(
    const Eigen::Matrix<var, R1, C1> &A,
    const Eigen::Matrix<double, R2, C2> &b) {
  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());

  check_square("mdivide_left_spd", "A", A);
  check_multiplicable("mdivide_left_spd", "A", A, "b", b);

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  internal::mdivide_left_spd_vd_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_spd_vd_vari<R1, C1, R2, C2>(A, b);

  res.vi() = Eigen::Map<matrix_vi>(&baseVari->variRefC_[0], b.rows(), b.cols());
  return res;
}

template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<var, R1, C2> mdivide_left_spd(
    const Eigen::Matrix<double, R1, C1> &A,
    const Eigen::Matrix<var, R2, C2> &b) {
  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());

  check_square("mdivide_left_spd", "A", A);
  check_multiplicable("mdivide_left_spd", "A", A, "b", b);

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  internal::mdivide_left_spd_dv_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_spd_dv_vari<R1, C1, R2, C2>(A, b);

  res.vi() = Eigen::Map<matrix_vi>(&baseVari->variRefC_[0], b.rows(), b.cols());

  return res;
}

}  // namespace math
}  // namespace stan
#endif
