#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_SPD_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <int R1, int C1, int R2, int C2>
class mdivide_left_spd_alloc : public chainable_alloc {
 public:
  virtual ~mdivide_left_spd_alloc() {}

  Eigen::LLT<Eigen::Matrix<double, R1, C1>> llt_;
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
    check_pos_definite("mdivide_left_spd", "A", alloc_->llt_);
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
    check_pos_definite("mdivide_left_spd", "A", alloc_->llt_);
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
    check_pos_definite("mdivide_left_spd", "A", alloc_->llt_);
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

template <
    typename EigMat1, typename EigMat2,
    require_all_eigen_matrix_base_vt<is_var, EigMat1, EigMat2> * = nullptr>
inline Eigen::Matrix<var, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_left_spd(const EigMat1 &A, const EigMat2 &b) {
  constexpr int R1 = EigMat1::RowsAtCompileTime;
  constexpr int C1 = EigMat1::ColsAtCompileTime;
  constexpr int R2 = EigMat2::RowsAtCompileTime;
  constexpr int C2 = EigMat2::ColsAtCompileTime;
  static constexpr const char *function = "mdivide_left_spd";
  check_multiplicable(function, "A", A, "b", b);
  const auto &A_ref = to_ref(A);
  check_symmetric(function, "A", A_ref);
  check_not_nan(function, "A", A_ref);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  internal::mdivide_left_spd_vv_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_spd_vv_vari<R1, C1, R2, C2>(A_ref, b);

  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(&baseVari->variRefC_[0], b.rows(), b.cols());
  return res;
}

template <typename EigMat1, typename EigMat2,
          require_eigen_matrix_base_vt<is_var, EigMat1> * = nullptr,
          require_eigen_matrix_base_vt<std::is_arithmetic, EigMat2> * = nullptr>
inline Eigen::Matrix<var, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_left_spd(const EigMat1 &A, const EigMat2 &b) {
  constexpr int R1 = EigMat1::RowsAtCompileTime;
  constexpr int C1 = EigMat1::ColsAtCompileTime;
  constexpr int R2 = EigMat2::RowsAtCompileTime;
  constexpr int C2 = EigMat2::ColsAtCompileTime;
  static constexpr const char *function = "mdivide_left_spd";
  check_multiplicable(function, "A", A, "b", b);
  const auto &A_ref = to_ref(A);
  check_symmetric(function, "A", A_ref);
  check_not_nan(function, "A", A_ref);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  internal::mdivide_left_spd_vd_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_spd_vd_vari<R1, C1, R2, C2>(A_ref, b);

  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(&baseVari->variRefC_[0], b.rows(), b.cols());
  return res;
}

template <typename EigMat1, typename EigMat2,
          require_eigen_matrix_base_vt<std::is_arithmetic, EigMat1> * = nullptr,
          require_eigen_matrix_base_vt<is_var, EigMat2> * = nullptr>
inline Eigen::Matrix<var, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_left_spd(const EigMat1 &A, const EigMat2 &b) {
  constexpr int R1 = EigMat1::RowsAtCompileTime;
  constexpr int C1 = EigMat1::ColsAtCompileTime;
  constexpr int R2 = EigMat2::RowsAtCompileTime;
  constexpr int C2 = EigMat2::ColsAtCompileTime;
  static constexpr const char *function = "mdivide_left_spd";
  check_multiplicable(function, "A", A, "b", b);
  const auto &A_ref = to_ref(A);
  check_symmetric(function, "A", A_ref);
  check_not_nan(function, "A", A_ref);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  // NOTE: this is not a memory leak, this vari is used in the
  // expression graph to evaluate the adjoint, but is not needed
  // for the returned matrix.  Memory will be cleaned up with the
  // arena allocator.
  internal::mdivide_left_spd_dv_vari<R1, C1, R2, C2> *baseVari
      = new internal::mdivide_left_spd_dv_vari<R1, C1, R2, C2>(A_ref, b);

  Eigen::Matrix<var, R1, C2> res(b.rows(), b.cols());
  res.vi() = Eigen::Map<matrix_vi>(&baseVari->variRefC_[0], b.rows(), b.cols());

  return res;
}

/**
 * Returns the solution of the system Ax=B where A is symmetric positive
 * definite.
 *
 * This overload handles arguments where one of T1 or T2 are
 * `var_value<T>` where `T` is an Eigen type. The other type can
 * also be a `var_value` or it can be a matrix type that inherits
 * from EigenBase
 *
 * @tparam T1 type of the first matrix
 * @tparam T2 type of the right-hand side matrix or vector
 *
 * @param A Matrix.
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if A is not square or B does not have
 * as many rows as A has columns.
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
          require_any_var_matrix_t<T1, T2> * = nullptr>
inline auto mdivide_left_spd(const T1 &A, const T2 &B) {
  using ret_val_type = plain_type_t<decltype(value_of(A) * value_of(B))>;
  using ret_type = var_value<ret_val_type>;

  if (A.size() == 0) {
    return ret_type(ret_val_type(0, B.cols()));
  }

  check_multiplicable("mdivide_left_spd", "A", A, "B", B);

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;

    check_symmetric("mdivide_left_spd", "A", arena_A.val());
    check_not_nan("mdivide_left_spd", "A", arena_A.val());

    auto A_llt = arena_A.val().llt();

    check_pos_definite("mdivide_left_spd", "A", A_llt);

    arena_t<Eigen::MatrixXd> arena_A_llt = A_llt.matrixL();
    arena_t<ret_type> res = A_llt.solve(arena_B.val());

    reverse_pass_callback([arena_A, arena_B, arena_A_llt, res]() mutable {
      promote_scalar_t<double, T2> adjB = res.adj();

      arena_A_llt.template triangularView<Eigen::Lower>().solveInPlace(adjB);
      arena_A_llt.template triangularView<Eigen::Lower>()
          .transpose()
          .solveInPlace(adjB);

      arena_A.adj() -= adjB * res.val_op().transpose();
      arena_B.adj() += adjB;
    });

    return ret_type(res);
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;

    check_symmetric("mdivide_left_spd", "A", arena_A.val());
    check_not_nan("mdivide_left_spd", "A", arena_A.val());

    auto A_llt = arena_A.val().llt();

    check_pos_definite("mdivide_left_spd", "A", A_llt);

    arena_t<Eigen::MatrixXd> arena_A_llt = A_llt.matrixL();
    arena_t<ret_type> res = A_llt.solve(value_of(B));

    reverse_pass_callback([arena_A, arena_A_llt, res]() mutable {
      promote_scalar_t<double, T2> adjB = res.adj();

      arena_A_llt.template triangularView<Eigen::Lower>().solveInPlace(adjB);
      arena_A_llt.template triangularView<Eigen::Lower>()
          .transpose()
          .solveInPlace(adjB);

      arena_A.adj() -= adjB * res.val().transpose().eval();
    });

    return ret_type(res);
  } else {
    const auto &A_ref = to_ref(value_of(A));
    arena_t<promote_scalar_t<var, T2>> arena_B = B;

    check_symmetric("mdivide_left_spd", "A", A_ref);
    check_not_nan("mdivide_left_spd", "A", A_ref);

    auto A_llt = A_ref.llt();

    check_pos_definite("mdivide_left_spd", "A", A_llt);

    arena_t<Eigen::MatrixXd> arena_A_llt = A_llt.matrixL();
    arena_t<ret_type> res = A_llt.solve(arena_B.val());

    reverse_pass_callback([arena_B, arena_A_llt, res]() mutable {
      promote_scalar_t<double, T2> adjB = res.adj();

      arena_A_llt.template triangularView<Eigen::Lower>().solveInPlace(adjB);
      arena_A_llt.template triangularView<Eigen::Lower>()
          .transpose()
          .solveInPlace(adjB);

      arena_B.adj() += adjB;
    });

    return ret_type(res);
  }
}

}  // namespace math
}  // namespace stan
#endif
