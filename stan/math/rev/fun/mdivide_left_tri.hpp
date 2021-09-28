#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_TRI_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_TRI_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

namespace internal {

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_vv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefA_;
  vari **variRefB_;
  vari **variRefC_;

  mdivide_left_tri_vv_vari(const Eigen::Matrix<var, R1, C1> &A,
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
                                                       * (A.rows() + 1) / 2))),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;

    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = j; i < M_; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    }

    Map<matrix_d> c_map(C_, M_, N_);
    Map<matrix_d> a_map(A_, M_, M_);
    a_map = A.val();
    c_map = B.val();
    Map<matrix_vi>(variRefB_, M_, N_) = B.vi();

    c_map = a_map.template triangularView<TriView>().solve(c_map);

    Map<matrix_vi>(variRefC_, M_, N_)
        = c_map.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjA;
    matrix_d adjB;

    adjB = Map<matrix_d>(A_, M_, M_)
               .template triangularView<TriView>()
               .transpose()
               .solve(Map<matrix_vi>(variRefC_, M_, N_).adj());
    adjA = -adjB * Map<matrix_d>(C_, M_, N_).transpose();

    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = j; i < adjA.rows(); i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    }
    Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_dv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefB_;
  vari **variRefC_;

  mdivide_left_tri_dv_vari(const Eigen::Matrix<double, R1, C1> &A,
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

    Map<matrix_d>(A_, M_, M_) = A;
    Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    Map<matrix_d> c_map(C_, M_, N_);
    c_map = B.val();

    c_map = Map<matrix_d>(A_, M_, M_)
                .template triangularView<TriView>()
                .solve(c_map);

    Map<matrix_vi>(variRefC_, M_, N_)
        = c_map.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;

    Map<matrix_vi>(variRefB_, M_, N_).adj()
        += Map<matrix_d>(A_, M_, M_)
               .template triangularView<TriView>()
               .transpose()
               .solve(Map<matrix_vi>(variRefC_, M_, N_).adj());
  }
};

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_vd_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefA_;
  vari **variRefC_;

  mdivide_left_tri_vd_vari(const Eigen::Matrix<var, R1, C1> &A,
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
                                                       * (A.rows() + 1) / 2))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;
    using Eigen::Matrix;

    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = j; i < M_; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    }
    Map<matrix_d> Ad(A_, M_, M_);
    Map<matrix_d> Cd(C_, M_, N_);
    Ad = A.val();

    Cd = Ad.template triangularView<TriView>().solve(B);

    Map<matrix_vi>(variRefC_, M_, N_)
        = Cd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    using Eigen::Matrix;
    Matrix<double, R1, C1> adjA(M_, M_);
    Matrix<double, R1, C2> adjC(M_, N_);

    adjC = Map<matrix_vi>(variRefC_, M_, N_).adj();

    adjA.noalias()
        = -Map<Matrix<double, R1, C1>>(A_, M_, M_)
               .template triangularView<TriView>()
               .transpose()
               .solve(adjC
                      * Map<Matrix<double, R1, C2>>(C_, M_, N_).transpose());

    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = j; i < adjA.rows(); i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    }
  }
};
}  // namespace internal

template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_eigen_vt<is_var, T1, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b) {
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  if (A.rows() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  auto *baseVari = new internal::mdivide_left_tri_vv_vari<
      TriView, T1::RowsAtCompileTime, T1::ColsAtCompileTime,
      T2::RowsAtCompileTime, T2::ColsAtCompileTime>(A, b);

  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> res(
      b.rows(), b.cols());
  res.vi()
      = Eigen::Map<matrix_vi>(&(baseVari->variRefC_[0]), b.rows(), b.cols());

  return res;
}
template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_eigen_vt<std::is_arithmetic, T1> * = nullptr,
          require_eigen_vt<is_var, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b) {
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  if (A.rows() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  auto *baseVari = new internal::mdivide_left_tri_dv_vari<
      TriView, T1::RowsAtCompileTime, T1::ColsAtCompileTime,
      T2::RowsAtCompileTime, T2::ColsAtCompileTime>(A, b);

  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> res(
      b.rows(), b.cols());
  res.vi()
      = Eigen::Map<matrix_vi>(&(baseVari->variRefC_[0]), b.rows(), b.cols());

  return res;
}
template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_eigen_vt<is_var, T1> * = nullptr,
          require_eigen_vt<std::is_arithmetic, T2> * = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b) {
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  if (A.rows() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  auto *baseVari = new internal::mdivide_left_tri_vd_vari<
      TriView, T1::RowsAtCompileTime, T1::ColsAtCompileTime,
      T2::RowsAtCompileTime, T2::ColsAtCompileTime>(A, b);

  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> res(
      b.rows(), b.cols());
  res.vi()
      = Eigen::Map<matrix_vi>(&(baseVari->variRefC_[0]), b.rows(), b.cols());

  return res;
}

/**
 * Returns the solution of the system Ax=B when A is triangular.
 *
 * This overload handles arguments where one of T1 or T2 are
 * `var_value<T>` where `T` is an Eigen type. The other type can
 * also be a `var_value` or it can be a matrix type that inherits
 * from EigenBase
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam T1 type of the triangular matrix
 * @tparam T2 type of the right-hand side matrix or vector
 *
 * @param A Triangular matrix.
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if A is not square or B does not have
 * as many rows as A has columns.
 */
template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_matrix_t<T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto mdivide_left_tri(const T1 &A, const T2 &B) {
  using ret_val_type = plain_type_t<decltype(value_of(A) * value_of(B))>;
  using ret_type = var_value<ret_val_type>;

  if (A.size() == 0) {
    return ret_type(ret_val_type(0, B.cols()));
  }

  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "B", B);

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    auto arena_A_val = to_arena(arena_A.val());

    arena_t<ret_type> res
        = arena_A_val.template triangularView<TriView>().solve(arena_B.val());

    reverse_pass_callback([arena_A, arena_B, arena_A_val, res]() mutable {
      promote_scalar_t<double, T2> adjB
          = arena_A_val.template triangularView<TriView>().transpose().solve(
              res.adj());

      arena_B.adj() += adjB;
      arena_A.adj() -= (adjB * res.val().transpose().eval())
                           .template triangularView<TriView>();
    });

    return ret_type(res);
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    auto arena_A_val = to_arena(arena_A.val());

    arena_t<ret_type> res
        = arena_A_val.template triangularView<TriView>().solve(value_of(B));

    reverse_pass_callback([arena_A, arena_A_val, res]() mutable {
      promote_scalar_t<double, T2> adjB
          = arena_A_val.template triangularView<TriView>().transpose().solve(
              res.adj());

      arena_A.adj() -= (adjB * res.val().transpose().eval())
                           .template triangularView<TriView>();
    });

    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<double, T1>> arena_A = value_of(A);
    arena_t<promote_scalar_t<var, T2>> arena_B = B;

    arena_t<ret_type> res
        = arena_A.template triangularView<TriView>().solve(arena_B.val());

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      promote_scalar_t<double, T2> adjB
          = arena_A.template triangularView<TriView>().transpose().solve(
              res.adj());

      arena_B.adj() += adjB;
    });

    return ret_type(res);
  }
}

}  // namespace math
}  // namespace stan
#endif
