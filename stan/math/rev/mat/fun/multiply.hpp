#ifndef STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/meta.hpp>
#include <boost/math/tools/promotion.hpp>
#include <type_traits>

namespace stan {
namespace math {
template <typename T1, typename T2, typename = void, typename = void,
          typename = void>
class multiply_mat_vari : public vari {};

/**
 * This is a subclass of the vari class for matrix
 * multiplication A * B where A is N by M and B
 * is M by K.
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if A or B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of A * B.
 *
 * @tparam Ta Scalar type for matrix A
 * @tparam Ra Rows for matrix A
 * @tparam Ca Columns for matrix A, Rows for matrix B
 * @tparam Tb Scalar type for matrix B
 * @tparam Cb Columns for matrix B
 */
template <typename T1, typename T2>
class multiply_mat_vari<T1, T2, enable_if_all_contains_var<T1, T2>,
                        enable_if_not_dot_product<T1, T2>, void> : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int B_cols_;
  int A_size_;
  int B_size_;
  double* Ad_;
  double* Bd_;
  vari** variRefA_;
  vari** variRefB_;
  vari** variRefAB_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A matrix
   * @param B matrix
   */
  multiply_mat_vari(const T1& A, const T2& B)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::instance_->memalloc_.alloc_array<double>(B_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefB_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(B_size_)),
        variRefAB_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A_rows_ * B_cols_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_vi>(variRefB_, A_cols_, B_cols_) = B.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_d> Bd(Bd_, A_cols_, B_cols_);
    Ad = A.val();
    Bd = B.val();

    Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
        = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjAB(A_rows_, B_cols_);

    adjAB = Map<matrix_vi>(variRefAB_, A_rows_, B_cols_).adj();
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
        += adjAB * Map<matrix_d>(Bd_, A_cols_, B_cols_).transpose();
    Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj()
        += Map<matrix_d>(Ad_, A_rows_, A_cols_).transpose() * adjAB;
  }
};

/**
 * This is a subclass of the vari class for matrix
 * multiplication A * B where A is 1 by M and B
 * is M by 1.
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if A or B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of A * B.
 *
 * @tparam Ta Scalar type for matrix A
 * @tparam Ca Columns for matrix A, Rows for matrix B
 * @tparam Tb Scalar type for matrix B
 */  // <Ta, 1, Ca, Tb, 1>
template <typename T1, typename T2>
class multiply_mat_vari<T1, T2, enable_if_all_contains_var<T1, T2>,
                        enable_if_dot_product<T1, T2>, void> : public vari {
 public:
  int size_;
  double* Ad_;
  double* Bd_;
  vari** variRefA_;
  vari** variRefB_;
  vari* variRefAB_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const T1& A, const T2& B)
      : vari(0.0),
        size_(A.cols()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(size_)),
        Bd_(ChainableStack::instance_->memalloc_.alloc_array<double>(size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(size_)),
        variRefB_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(size_)) {
    using Eigen::Map;
    Map<row_vector_vi>(variRefA_, size_) = A.vi();
    Map<vector_vi>(variRefB_, size_) = B.vi();
    Map<row_vector_d> Ad(Ad_, size_);
    Map<vector_d> Bd(Bd_, size_);
    Ad = A.val();
    Bd = B.val();

    variRefAB_ = new vari(Ad * Bd, false);
  }

  virtual void chain() {
    using Eigen::Map;

    double adjAB = variRefAB_->adj_;
    Map<vector_vi>(variRefA_, size_).adj() += adjAB * Map<vector_d>(Bd_, size_);
    Map<vector_vi>(variRefB_, size_).adj() += Map<vector_d>(Ad_, size_) * adjAB;
  }
};

/**
 * This is a subclass of the vari class for matrix
 * multiplication A * B where A is an N by M
 * matrix of double and B is M by K.
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if A or B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of A * B.
 *
 * @tparam Ra Rows for matrix A
 * @tparam Ca Columns for matrix A, Rows for matrix B
 * @tparam Tb Scalar type for matrix B
 * @tparam Cb Columns for matrix B
 */  // <double, Ra, Ca, Tb, Cb>
template <typename T1, typename T2>
class multiply_mat_vari<T1, T2, enable_if_contains_arithmetic<T1>,
                        enable_if_contains_var<T2>,
                        enable_if_not_dot_product<T1, T2>> : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int B_cols_;
  int A_size_;
  int B_size_;
  double* Ad_;
  double* Bd_;
  vari** variRefB_;
  vari** variRefAB_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const T1& A, const T2& B)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::instance_->memalloc_.alloc_array<double>(B_size_)),
        variRefB_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(B_size_)),
        variRefAB_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A_rows_ * B_cols_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefB_, A_cols_, B_cols_) = B.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_d> Bd(Bd_, A_cols_, B_cols_);
    Ad = A;
    Bd = B.val();

    Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
        = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjAB = Map<matrix_vi>(variRefAB_, A_rows_, B_cols_).adj();

    Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj()
        += Map<matrix_d>(Ad_, A_rows_, A_cols_).transpose() * adjAB;
  }
};

/**
 * This is a subclass of the vari class for matrix
 * multiplication A * B where A is a double
 * row vector of length M and B is a vector of
 * length M.
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if A or B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of A * B.
 *
 * @tparam Ca Columns for matrix A, Rows for matrix B
 * @tparam Tb Scalar type for matrix B
 */  // <double, 1, Ca, Tb, 1>
template <typename T1, typename T2>
class multiply_mat_vari<T1, T2, enable_if_contains_arithmetic<T1>,
                        enable_if_contains_var<T2>,
                        enable_if_dot_product<T1, T2>> : public vari {
 public:
  int size_;
  double* Ad_;
  double* Bd_;
  vari** variRefB_;
  vari* variRefAB_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const T1& A, const T2& B)
      : vari(0.0),
        size_(A.cols()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(size_)),
        Bd_(ChainableStack::instance_->memalloc_.alloc_array<double>(size_)),
        variRefB_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(size_)) {
    using Eigen::Map;
    Map<row_vector_d> Ad(Ad_, size_);
    Map<vector_d> Bd(Bd_, size_);
    Map<vector_vi>(variRefB_, size_) = B.vi();
    Ad = A;
    Bd = B.val();
    variRefAB_ = new vari(Ad * Bd, false);
  }

  virtual void chain() {
    using Eigen::Map;
    Map<vector_vi>(variRefB_, size_).adj()
        += Map<vector_d>(Ad_, size_) * variRefAB_->adj_;
  }
};

/**
 * This is a subclass of the vari class for matrix
 * multiplication A * B where A is N by M and B
 * is an M by K matrix of doubles.
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if A or B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of A * B.
 *
 * @tparam Ta Scalar type for matrix A
 * @tparam Ra Rows for matrix A
 * @tparam Ca Columns for matrix A, Rows for matrix B
 * @tparam Cb Columns for matrix B
 */  // <Ta, Ra, Ca, double, Cb>
template <typename T1, typename T2>
class multiply_mat_vari<T1, T2, enable_if_contains_var<T1>,
                        enable_if_contains_arithmetic<T2>,
                        enable_if_not_dot_product<T1, T2>> : public vari {
 public:
  int A_rows_;
  int A_cols_;
  int B_cols_;
  int A_size_;
  int B_size_;
  double* Ad_;
  double* Bd_;
  vari** variRefA_;
  vari** variRefAB_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const T1& A, const T2& B)
      : vari(0.0),
        A_rows_(A.rows()),
        A_cols_(A.cols()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::instance_->memalloc_.alloc_array<double>(B_size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(A_size_)),
        variRefAB_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A_rows_ * B_cols_)) {
    using Eigen::Map;
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_) = A.vi();
    Map<matrix_d> Ad(Ad_, A_rows_, A_cols_);
    Map<matrix_d> Bd(Bd_, A_cols_, B_cols_);
    Ad = A.val();
    Bd = B.val();

    Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
        = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjAB = Map<matrix_vi>(variRefAB_, A_rows_, B_cols_).adj();

    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
        += adjAB * Map<matrix_d>(Bd_, A_cols_, B_cols_).transpose();
  }
};

/**
 * This is a subclass of the vari class for matrix
 * multiplication A * B where A is a row
 * vector of length M and B is a vector of length M
 * of doubles.
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if A or B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of A * B.
 *
 * @tparam Ta Scalar type for matrix A
 * @tparam Ra Rows for matrix A
 * @tparam Ca Columns for matrix A, Rows for matrix B
 * @tparam Tb Scalar type for matrix B
 * @tparam Cb Columns for matrix B
 */  // <Ta, 1, Ca, double, 1>
template <typename T1, typename T2>
class multiply_mat_vari<T1, T2, enable_if_contains_var<T1>,
                        enable_if_contains_arithmetic<T2>,
                        enable_if_dot_product<T1, T2>> : public vari {
 public:
  int size_;
  double* Ad_;
  double* Bd_;
  vari** variRefA_;
  vari* variRefAB_;

  /**
   * Constructor for multiply_mat_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const T1& A, const T2& B)
      : vari(0.0),
        size_(A.cols()),
        Ad_(ChainableStack::instance_->memalloc_.alloc_array<double>(size_)),
        Bd_(ChainableStack::instance_->memalloc_.alloc_array<double>(size_)),
        variRefA_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(size_)) {
    using Eigen::Map;
    Map<row_vector_vi>(variRefA_, size_) = A.vi();
    Map<row_vector_d> Ad(Ad_, size_);
    Map<vector_d> Bd(Bd_, size_);
    Ad = A.val();
    Bd = B;

    variRefAB_ = new vari(Ad * Bd, false);
  }

  virtual void chain() {
    using Eigen::Map;

    Map<row_vector_vi>(variRefA_, size_).adj()
        += variRefAB_->adj_ * Map<row_vector_d>(Bd_, size_);
  }
};

/**
 * Return the product of scalar and matrix.
 * @tparam T1 Type of LHS input, either scalar or matrix.
 * @tparam T2 Type of RHS input, either scalar or matrix.
 * @param[in] A Specified scalar
 * @param[in] B Matrix
 * @return Product of scalar and matrix
 */
template <typename T1, typename T2,
          enable_if_all_eigen_or_stan_scalar<T1, T2>* = nullptr,
          enable_if_all_not_eigen<T1, T2>* = nullptr,
          enable_if_any_contains_var<T1, T2>* = nullptr,
          enable_if_not_dot_product<T1, T2>* = nullptr>
inline auto multiply(const T1& A, const T2& B) {
  return A * B;
}

/**
 * Return the product of two matrices. Specialized for eigen matrices of vars that are not dot products
 * @tparam T1 Type of LHS Eigen Matrix
 * @tparam T2 Type of RHS Eigen Matrix
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return Product of scalar and matrix.
 */
template <typename T1, typename T2, enable_if_all_eigen<T1, T2>* = nullptr,
          enable_if_any_contains_var<T1, T2>* = nullptr,
          enable_if_not_dot_product<T1, T2>* = nullptr>
inline auto multiply(const T1& A, const T2& B) {
  check_multiplicable("multiply", "A", A, "B", B);
  check_not_nan("multiply", "A", A);
  check_not_nan("multiply", "B", B);

  // Memory managed with the arena allocator.
  multiply_mat_vari<T1, T2>* baseVari = new multiply_mat_vari<T1, T2>(A, B);
  Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime> AB_v(
      A.rows(), B.cols());
  AB_v.vi()
      = Eigen::Map<matrix_vi>(&baseVari->variRefAB_[0], A.rows(), B.cols());

  return AB_v;
}

/**
 * Return the scalar product of a row vector and column vector.
 * @tparam T1 Type of Eigen Matrix
 * @tparam T2 Type of Eigen Matrix
 * @param[in] A Row vector
 * @param[in] B Column vector
 * @return Scalar product of row vector and vector
 */
template <typename T1, typename T2, enable_if_all_eigen<T1, T2>* = nullptr,
          enable_if_any_contains_var<T1, T2>* = nullptr,
          enable_if_dot_product<T1, T2>* = nullptr>
inline auto multiply(const T1& A, const T2& B) {
  check_multiplicable("multiply", "A", A, "B", B);
  check_not_nan("multiply", "A", A);
  check_not_nan("multiply", "B", B);
  // Memory managed with the arena allocator.
  multiply_mat_vari<T1, T2>* baseVari = new multiply_mat_vari<T1, T2>(A, B);
  var AB_v;

  AB_v.vi_ = baseVari->variRefAB_;
  return AB_v;
}

}  // namespace math
}  // namespace stan
#endif
