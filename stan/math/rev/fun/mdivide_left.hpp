#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <int R1, int C1, int R2, int C2>
class mdivide_left_vv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefA_;
  vari **variRefB_;
  vari **variRefC_;

  mdivide_left_vv_vari(const Eigen::Matrix<var, R1, C1> &A,
                       const Eigen::Matrix<var, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        C_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * B.rows()
                                                       * B.cols()))),
        variRefA_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;

    Map<matrix_d> Ad(A_, M_, M_);
    Map<matrix_d> Cd(C_, M_, N_);
    Ad = A.val();
    Cd = Ad.colPivHouseholderQr().solve(B.val());

    Map<matrix_vi>(variRefA_, M_, M_) = A.vi();
    Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    Map<matrix_vi>(variRefC_, M_, N_)
        = Cd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjB = Map<matrix_d>(A_, M_, M_)
                        .transpose()
                        .colPivHouseholderQr()
                        .solve(Map<matrix_vi>(variRefC_, M_, N_).adj());

    Map<matrix_vi>(variRefA_, M_, M_).adj()
        -= adjB * Map<matrix_d>(C_, M_, N_).transpose();
    Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

template <int R1, int C1, int R2, int C2>
class mdivide_left_dv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefB_;
  vari **variRefC_;

  mdivide_left_dv_vari(const Eigen::Matrix<double, R1, C1> &A,
                       const Eigen::Matrix<var, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        C_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * B.rows()
                                                       * B.cols()))),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;

    Map<matrix_d> Ad(A_, M_, M_);
    Map<matrix_d> Cd(C_, M_, N_);
    Ad = A;
    Cd = Ad.colPivHouseholderQr().solve(B.val());
    Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    Map<matrix_vi>(variRefC_, M_, N_)
        = Cd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;

    Map<matrix_vi>(variRefB_, M_, N_).adj()
        += Map<matrix_d>(A_, M_, M_)
               .transpose()
               .colPivHouseholderQr()
               .solve(Map<matrix_vi>(variRefC_, M_, N_).adj());
  }
};

template <int R1, int C1, int R2, int C2>
class mdivide_left_vd_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefA_;
  vari **variRefC_;

  mdivide_left_vd_vari(const Eigen::Matrix<var, R1, C1> &A,
                       const Eigen::Matrix<double, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        C_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * B.rows()
                                                       * B.cols()))),
        variRefA_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * A.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;

    Map<matrix_vi>(variRefA_, M_, M_) = A.vi();
    Map<matrix_d> Ad(A_, M_, M_);
    Map<matrix_d> Cd(C_, M_, N_);

    Ad = A.val();
    Cd = Ad.colPivHouseholderQr().solve(B);
    Map<matrix_vi>(variRefC_, M_, N_)
        = Cd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;

    matrix_d adjC = Map<matrix_vi>(variRefC_, M_, N_).adj();

    Map<matrix_vi>(variRefA_, M_, M_).adj()
        -= Map<matrix_d>(A_, M_, M_)
               .transpose()
               .colPivHouseholderQr()
               .solve(adjC * Map<matrix_d>(C_, M_, N_).transpose());
  }
};
}  // namespace internal

template <typename T1, typename T2,
          require_all_eigen_vt<is_var, T1, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left(const T1 &A, const T2 &b) {
  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  auto *baseVari = new internal::mdivide_left_vv_vari<
      T1::RowsAtCompileTime, T1::ColsAtCompileTime, T2::RowsAtCompileTime,
      T2::ColsAtCompileTime>(A, b);

  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> res(
      b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

template <typename T1, typename T2, require_eigen_vt<is_var, T1> * = nullptr,
          require_eigen_vt<std::is_arithmetic, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left(const T1 &A, const T2 &b) {
  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  auto *baseVari = new internal::mdivide_left_vd_vari<
      T1::RowsAtCompileTime, T1::ColsAtCompileTime, T2::RowsAtCompileTime,
      T2::ColsAtCompileTime>(A, b);

  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> res(
      b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

template <typename T1, typename T2,
          require_eigen_vt<std::is_arithmetic, T1> * = nullptr,
          require_eigen_vt<is_var, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left(const T1 &A, const T2 &b) {
  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  auto *baseVari = new internal::mdivide_left_dv_vari<
      T1::RowsAtCompileTime, T1::ColsAtCompileTime, T2::RowsAtCompileTime,
      T2::ColsAtCompileTime>(A, b);

  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> res(
      b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(baseVari->variRefC_, res.rows(), res.cols());

  return res;
}

}  // namespace math
}  // namespace stan
#endif
