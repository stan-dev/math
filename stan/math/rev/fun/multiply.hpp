#ifndef STAN_MATH_REV_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim.hpp>
#include <type_traits>

namespace stan {
namespace math {

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
 * @tparam Ta type of elements in matrix A
 * @tparam Ra number of rows in matrix A, can be Eigen::Dynamic
 * @tparam Ca number of columns in matrix A and rows in matrix B
 * @tparam Tb type of elements in matrix B
 * @tparam Cb number of columns in matrix B, can be Eigen::Dynamic
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Cb>
class multiply_mat_vari : public vari {
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
   * controlled by the second argument to
   * vari's constructor.
   *
   * @param A matrix
   * @param B matrix
   */
  multiply_mat_vari(const Eigen::Matrix<Ta, Ra, Ca>& A,
                    const Eigen::Matrix<Tb, Ca, Cb>& B)
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
#ifdef STAN_OPENCL
    if (Ad.rows() * Ad.cols() * Bd.cols()
        > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
      matrix_cl<double> Ad_cl(Ad);
      matrix_cl<double> Bd_cl(Bd);
      matrix_cl<double> variRefAB_cl = Ad_cl * Bd_cl;
      matrix_d temp = from_matrix_cl(variRefAB_cl);
      Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
          = temp.unaryExpr([](double x) { return new vari(x, false); });
    } else {
      Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
          = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
    }
#else
    Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
        = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
#endif
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjAB(A_rows_, B_cols_);
    adjAB = Map<matrix_vi>(variRefAB_, A_rows_, B_cols_).adj();
#ifdef STAN_OPENCL
    if (A_rows_ * A_cols_ * B_cols_
        > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
      matrix_cl<double> adjAB_cl(adjAB);
      matrix_cl<double> Ad_cl(Ad_, A_rows_, A_cols_);
      matrix_cl<double> Bd_cl(Bd_, A_cols_, B_cols_);
      matrix_cl<double> variRefA_cl = adjAB_cl * transpose(Bd_cl);
      matrix_cl<double> variRefB_cl = transpose(Ad_cl) * adjAB_cl;
      matrix_d temp_variRefA = from_matrix_cl(variRefA_cl);
      matrix_d temp_variRefB = from_matrix_cl(variRefB_cl);
      Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj() += temp_variRefA;
      Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj() += temp_variRefB;
    } else {
      Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
          += adjAB * Map<matrix_d>(Bd_, A_cols_, B_cols_).transpose();
      Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj()
          += Map<matrix_d>(Ad_, A_rows_, A_cols_).transpose() * adjAB;
    }
#else
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
        += adjAB * Map<matrix_d>(Bd_, A_cols_, B_cols_).transpose();
    Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj()
        += Map<matrix_d>(Ad_, A_rows_, A_cols_).transpose() * adjAB;
#endif
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
 * @tparam Ta type of elements in matrix A
 * @tparam Ca number of columns in matrix A and rows in matrix B
 * @tparam Tb type of elements in matrix B
 */
template <typename Ta, int Ca, typename Tb>
class multiply_mat_vari<Ta, 1, Ca, Tb, 1> : public vari {
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
   * controlled by the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const Eigen::Matrix<Ta, 1, Ca>& A,
                    const Eigen::Matrix<Tb, Ca, 1>& B)
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
 * @tparam Ra number of rows in matrix A, can be Eigen::Dynamic
 * @tparam Ca number of columns in matrix A and rows in matrix B
 * @tparam Tb type of elements in matrix B
 * @tparam Cb number of columns in matrix B, can be Eigen::Dynamic
 */
template <int Ra, int Ca, typename Tb, int Cb>
class multiply_mat_vari<double, Ra, Ca, Tb, Cb> : public vari {
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
   * controlled by the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const Eigen::Matrix<double, Ra, Ca>& A,
                    const Eigen::Matrix<Tb, Ca, Cb>& B)
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
#ifdef STAN_OPENCL
    if (Ad.rows() * Ad.cols() * Bd.cols()
        > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
      matrix_cl<double> Ad_cl(Ad);
      matrix_cl<double> Bd_cl(Bd);
      matrix_cl<double> variRefAB_cl = Ad_cl * Bd_cl;
      matrix_d temp = from_matrix_cl(variRefAB_cl);
      Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
          = temp.unaryExpr([](double x) { return new vari(x, false); });
    } else {
      Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
          = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
    }
#else
    Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
        = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
#endif
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjAB = Map<matrix_vi>(variRefAB_, A_rows_, B_cols_).adj();
#ifdef STAN_OPENCL
    if (A_rows_ * A_cols_ * B_cols_
        > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
      matrix_cl<double> adjAB_cl(adjAB);
      matrix_cl<double> Ad_cl(Ad_, A_rows_, A_cols_);
      matrix_cl<double> variRefB_cl = transpose(Ad_cl) * adjAB_cl;
      matrix_d temp_variRefB = from_matrix_cl(variRefB_cl);
      Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj() += temp_variRefB;
    } else {
      Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj()
          += Map<matrix_d>(Ad_, A_rows_, A_cols_).transpose() * adjAB;
    }
#else
    Map<matrix_vi>(variRefB_, A_cols_, B_cols_).adj()
        += Map<matrix_d>(Ad_, A_rows_, A_cols_).transpose() * adjAB;
#endif
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
 * @tparam Ca number of columns in matrix A and rows in matrix B
 * @tparam Tb type of elements in matrix B
 */
template <int Ca, typename Tb>
class multiply_mat_vari<double, 1, Ca, Tb, 1> : public vari {
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
   * controlled by the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const Eigen::Matrix<double, 1, Ca>& A,
                    const Eigen::Matrix<Tb, Ca, 1>& B)
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
 * @tparam Ta type of elements in matrix A
 * @tparam Ra number of rows in matrix A, can be Eigen::Dynamic
 * @tparam Ca number of columns in matrix A and rows in matrix B
 * @tparam Cb number of columns in matrix B, can be Eigen::Dynamic
 */
template <typename Ta, int Ra, int Ca, int Cb>
class multiply_mat_vari<Ta, Ra, Ca, double, Cb> : public vari {
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
   * controlled by the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const Eigen::Matrix<Ta, Ra, Ca>& A,
                    const Eigen::Matrix<double, Ca, Cb>& B)
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
#ifdef STAN_OPENCL
    if (Ad.rows() * Ad.cols() * Bd.cols()
        > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
      matrix_cl<double> Ad_cl(Ad);
      matrix_cl<double> Bd_cl(Bd);
      matrix_cl<double> variRefAB_cl = Ad_cl * Bd_cl;
      matrix_d temp = from_matrix_cl(variRefAB_cl);
      Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
          = temp.unaryExpr([](double x) { return new vari(x, false); });
    } else {
      Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
          = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
    }
#else
    Map<matrix_vi>(variRefAB_, A_rows_, B_cols_)
        = (Ad * Bd).unaryExpr([](double x) { return new vari(x, false); });
#endif
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjAB = Map<matrix_vi>(variRefAB_, A_rows_, B_cols_).adj();
#ifdef STAN_OPENCL
    if (A_rows_ * A_cols_ * B_cols_
        > opencl_context.tuning_opts().multiply_dim_prod_worth_transfer) {
      matrix_cl<double> adjAB_cl(adjAB);
      matrix_cl<double> Bd_cl(Bd_, A_cols_, B_cols_);
      matrix_cl<double> variRefA_cl = adjAB_cl * transpose(Bd_cl);
      matrix_d temp_variRefA = from_matrix_cl(variRefA_cl);
      Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj() += temp_variRefA;
    } else {
      Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
          += adjAB * Map<matrix_d>(Bd_, A_cols_, B_cols_).transpose();
    }
#else
    Map<matrix_vi>(variRefA_, A_rows_, A_cols_).adj()
        += adjAB * Map<matrix_d>(Bd_, A_cols_, B_cols_).transpose();
#endif
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
 * @tparam Ta type of elements in matrix A
 * @tparam Ca number of columns in matrix A and rows in matrix B
 */
template <typename Ta, int Ca>
class multiply_mat_vari<Ta, 1, Ca, double, 1> : public vari {
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
   * controlled by the second argument to
   * vari's constructor.
   *
   * @param A row vector
   * @param B vector
   */
  multiply_mat_vari(const Eigen::Matrix<Ta, 1, Ca>& A,
                    const Eigen::Matrix<double, Ca, 1>& B)
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
 * Return the product of two matrices.
 *
 * @tparam Ta type of elements in matrix A
 * @tparam Ra number of rows in matrix A, can be Eigen::Dynamic
 * @tparam Ca number of columns in matrix A and rows in matrix B
 * @tparam Tb type of elements in matrix B
 * @tparam Cb number of columns in matrix B, can be Eigen::Dynamic
 *
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return Product of scalar and matrix.
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_any_eigen_vt<is_var, Mat1, Mat2>* = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2>* = nullptr>
inline auto multiply(const Mat1& A, const Mat2& B) {
  using Ta = value_type_t<Mat1>;
  using Tb = value_type_t<Mat2>;
  constexpr int Ra = Mat1::RowsAtCompileTime;
  constexpr int Ca = Mat1::ColsAtCompileTime;
  constexpr int Cb = Mat2::ColsAtCompileTime;
  check_size_match("multiply", "Columns of A", A.cols(), "Rows of B", B.rows());
  check_not_nan("multiply", "A", A);
  check_not_nan("multiply", "B", B);

  // Memory managed with the arena allocator.
  multiply_mat_vari<Ta, Ra, Ca, Tb, Cb>* baseVari
      = new multiply_mat_vari<Ta, Ra, Ca, Tb, Cb>(A, B);
  Eigen::Matrix<var, Ra, Cb> AB_v(A.rows(), B.cols());
  AB_v.vi()
      = Eigen::Map<matrix_vi>(&baseVari->variRefAB_[0], A.rows(), B.cols());

  return AB_v;
}

/**
 * Return the scalar product of a row vector and
 * a vector.
 *
 * @tparam RowVec type of row vector A
 * @tparam ColVec type of column vector B
 *
 * @param[in] A Row vector
 * @param[in] B Column vector
 * @return Scalar product of row vector and vector
 */
template <
    typename RowVec, typename ColVec,
    require_any_var_t<value_type_t<RowVec>, value_type_t<ColVec>>* = nullptr,
    require_eigen_row_and_col_t<RowVec, ColVec>* = nullptr>
inline var multiply(const RowVec& A, const ColVec& B) {
  using RowVecScalar = value_type_t<RowVec>;
  using ColVecScalar = value_type_t<ColVec>;
  constexpr int Ca = RowVec::ColsAtCompileTime;
  check_matching_sizes("multiply", "A", A, "B", B);
  check_not_nan("multiply", "A", A);
  check_not_nan("multiply", "B", B);

  // Memory managed with the arena allocator.
  multiply_mat_vari<RowVecScalar, 1, Ca, ColVecScalar, 1>* baseVari
      = new multiply_mat_vari<RowVecScalar, 1, Ca, ColVecScalar, 1>(A, B);
  var AB_v;
  AB_v.vi_ = baseVari->variRefAB_;
  return AB_v;
}

}  // namespace math
}  // namespace stan
#endif
