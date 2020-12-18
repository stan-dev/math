#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_alloc.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/core/typedefs.hpp>
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
  std::shared_ptr<Eigen::LDLT<Eigen::Matrix<double, R1, C1>>> ldltP_;
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
    const LDLT_factor<var, R1, C1> &A, const EigMat &b) {
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
    const LDLT_factor<var, R1, C1> &A, const EigMat &b) {
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
    const LDLT_factor<double, R1, C1> &A, const EigMat &b) {
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

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T type of B
 * @param A LDLT_factor2
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
template <typename T1, bool alloc_in_arena, typename T2,
          require_all_matrix_t<T1, T2> * = nullptr,
          require_any_st_var<T1, T2> * = nullptr>
inline auto mdivide_left_ldlt(LDLT_factor2<T1, alloc_in_arena> &A,
                              const T2 &B) {
  using ret_val_type
      = Eigen::Matrix<double, Eigen::Dynamic, T2::ColsAtCompileTime>;
  using ret_type = promote_var_matrix_t<ret_val_type, T1, T2>;

  check_multiplicable("mdivide_left_ldlt", "A", A.matrix().val(), "B", B);

  const auto &B_ref = to_ref(B);

  if (A.matrix().size() == 0) {
    return ret_type(ret_val_type(0, B.cols()));
  }

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    arena_t<ret_type> res = A.ldlt().solve(arena_B.val());

    reverse_pass_callback([A, arena_B, res]() mutable {
      promote_scalar_t<double, T2> adjB = A.ldlt().solve(res.adj());

      forward_as<promote_scalar_t<var, T1>>(A.matrix()).adj()
          -= adjB * res.val().transpose().eval();
      arena_B.adj() += adjB;
    });

    return ret_type(res);
  } else if (!is_constant<T1>::value) {
    arena_t<ret_type> res = A.ldlt().solve(value_of(B));

    reverse_pass_callback([A, res]() mutable {
      promote_scalar_t<double, T2> adjB = A.ldlt().solve(res.adj());

      forward_as<promote_scalar_t<var, T1>>(A.matrix()).adj()
          -= adjB * res.val().transpose().eval();
    });

    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    arena_t<ret_type> res = A.ldlt().solve(arena_B.val());

    reverse_pass_callback([A, arena_B, res]() mutable {
      promote_scalar_t<double, T2> adjB = A.ldlt().solve(res.adj());

      arena_B.adj() += adjB;
    });

    return ret_type(res);
  }
}

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T type of B
 * @tparam R rowtype of A
 * @tparam C coltype of A
 * @param A LDLT_factor
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
/*template <typename T, int R, int C,
          require_var_matrix_t<T>* = nullptr>
inline auto mdivide_left_ldlt(const LDLT_factor<double, R, C>& A, const T& B) {
  using ret_val_type = Eigen::Matrix<double, R, T::ColsAtCompileTime>;
  using ret_type = var_value<ret_val_type>;

  check_multiplicable("mdivide_left_ldlt", "A", A, "B", B.val());

  if (A.rows() == 0) {
    return ret_type(ret_val_type(0, B.cols()));
  }

  arena_t<ret_type> res = A.solve(B.val());

  reverse_pass_callback([A, B, res]() mutable {
    B.adj() += A.solve(res.adj());
  });

  return ret_type(res);
  }*/

}  // namespace math
}  // namespace stan
#endif
