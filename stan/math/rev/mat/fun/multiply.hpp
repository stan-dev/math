#ifndef STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

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
     * @tparam TA Scalar type for matrix A
     * @tparam RA Rows for matrix A
     * @tparam CA Columns for matrix A, Rows for matrix B
     * @tparam TB Scalar type for matrix B
     * @tparam CB Columns for matrix B
     */
    template <typename TA, int RA, int CA, typename TB, int CB>
    class multiply_mat_vari : public vari {
      public:
        int A_rows_, A_cols_, B_cols_, A_size_, B_size_;
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
        multiply_mat_vari(const Eigen::Matrix<TA, RA, CA>& A,
                          const Eigen::Matrix<TB, CA, CB>& B)
          : vari(0.0),
          A_rows_(A.rows()), A_cols_(A.cols()),
          B_cols_(B.cols()), A_size_(A.size()), B_size_(B.size()),
          Ad_(ChainableStack::memalloc_.alloc_array<double>(A_size_)),
          Bd_(ChainableStack::memalloc_.alloc_array<double>(B_size_)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>(A_size_)),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>(B_size_)),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>(A_rows_
                                                                  * B_cols_)) {
            using Eigen::Map;
            using Eigen::MatrixXd;
            for (size_type i = 0; i < A.size(); ++i) {
              variRefA_[i] = A.coeffRef(i).vi_;
              Ad_[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < B.size(); ++i) {
              variRefB_[i] = B.coeffRef(i).vi_;
              Bd_[i] = B.coeffRef(i).val();
            }
            MatrixXd AB
              = Map<MatrixXd>(Ad_, A_rows_, A_cols_)
              * Map<MatrixXd>(Bd_, A_cols_, B_cols_);
            for (size_type i = 0; i < AB.size(); ++i)
              variRefAB_[i] = new vari(AB.coeffRef(i), false);
          }

        virtual void chain() {
          using Eigen::MatrixXd;
          using Eigen::Map;
          MatrixXd adjAB(A_rows_, B_cols_);
          MatrixXd adjA(A_rows_, A_cols_);
          MatrixXd adjB(A_cols_, B_cols_);

          for (size_type i = 0; i < adjAB.size(); ++i)
            adjAB(i) = variRefAB_[i]->adj_;
          adjA = adjAB
            * Map<MatrixXd>(Bd_, A_cols_, B_cols_).transpose();
          adjB = Map<MatrixXd>(Ad_, A_rows_, A_cols_).transpose()
            * adjAB;
          for (size_type i = 0; i < A_size_; ++i)
            variRefA_[i]->adj_ += adjA(i);
          for (size_type i = 0; i < B_size_; ++i)
            variRefB_[i]->adj_ += adjB(i);
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
     * @tparam TA Scalar type for matrix A
     * @tparam CA Columns for matrix A, Rows for matrix B
     * @tparam TB Scalar type for matrix B
     */
    template <typename TA, int CA, typename TB>
    class multiply_mat_vari<TA, 1, CA, TB, 1> : public vari {
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
        multiply_mat_vari(const Eigen::Matrix<TA, 1, CA>& A,
                          const Eigen::Matrix<TB, CA, 1>& B)
          : vari(0.0),
          size_(A.cols()),
          Ad_(ChainableStack::memalloc_.alloc_array<double>(size_)),
          Bd_(ChainableStack::memalloc_.alloc_array<double>(size_)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>(size_)),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {
            using Eigen::Map;
            using Eigen::VectorXd;
            using Eigen::RowVectorXd;
            for (size_type i = 0; i < size_; ++i) {
              variRefA_[i] = A.coeffRef(i).vi_;
              Ad_[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < size_; ++i) {
              variRefB_[i] = B.coeffRef(i).vi_;
              Bd_[i] = B.coeffRef(i).val();
            }
            double AB = Map<RowVectorXd>(Ad_, 1, size_)
                        * Map<VectorXd>(Bd_, size_, 1);
            variRefAB_ = new vari(AB, false);
          }

        virtual void chain() {
          using Eigen::VectorXd;
          using Eigen::RowVectorXd;
          using Eigen::Map;
          double adjAB;
          RowVectorXd adjA(size_);
          VectorXd adjB(size_);

          adjAB = variRefAB_->adj_;
          adjA = adjAB
            * Map<VectorXd>(Bd_, size_, 1).transpose();
          adjB = Map<RowVectorXd>(Ad_, 1, size_).transpose()
            * adjAB;
          for (size_type i = 0; i < size_; ++i)
            variRefA_[i]->adj_ += adjA(i);
          for (size_type i = 0; i < size_; ++i)
            variRefB_[i]->adj_ += adjB(i);
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
     * @tparam RA Rows for matrix A
     * @tparam CA Columns for matrix A, Rows for matrix B
     * @tparam TB Scalar type for matrix B
     * @tparam CB Columns for matrix B
     */
    template <int RA, int CA, typename TB, int CB>
    class multiply_mat_vari<double, RA, CA, TB, CB> : public vari {
      public:
        int A_rows_, A_cols_, B_cols_, A_size_, B_size_;
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
        multiply_mat_vari(const Eigen::Matrix<double, RA, CA>& A,
                          const Eigen::Matrix<TB, CA, CB>& B)
          : vari(0.0),
          A_rows_(A.rows()), A_cols_(A.cols()),
          B_cols_(B.cols()), A_size_(A.size()), B_size_(B.size()),
          Ad_(ChainableStack::memalloc_.alloc_array<double>(A_size_)),
          Bd_(ChainableStack::memalloc_.alloc_array<double>(B_size_)),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>(B_size_)),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>(A_rows_
                                                                  * B_cols_)) {
            using Eigen::MatrixXd;
            using Eigen::Map;
            for (size_type i = 0; i < A.size(); ++i)
              Ad_[i] = A.coeffRef(i);
            for (size_type i = 0; i < B.size(); ++i) {
              variRefB_[i] = B.coeffRef(i).vi_;
              Bd_[i] = B.coeffRef(i).val();
            }
            MatrixXd AB
              = Map<MatrixXd>(Ad_, A_rows_, A_cols_)
              * Map<MatrixXd>(Bd_, A_cols_, B_cols_);
            for (size_type i = 0; i < AB.size(); ++i)
              variRefAB_[i] = new vari(AB.coeffRef(i), false);
          }

        virtual void chain() {
          using Eigen::MatrixXd;
          using Eigen::Map;
          MatrixXd adjAB(A_rows_, B_cols_);
          MatrixXd adjB(A_cols_, B_cols_);

          for (size_type i = 0; i < adjAB.size(); ++i)
            adjAB(i) = variRefAB_[i]->adj_;
          adjB = Map<MatrixXd>(Ad_, A_rows_, A_cols_).transpose()
            * adjAB;
          for (size_type i = 0; i < B_size_; ++i)
            variRefB_[i]->adj_ += adjB(i);
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
     * @tparam CA Columns for matrix A, Rows for matrix B
     * @tparam TB Scalar type for matrix B
     */
    template <int CA, typename TB>
    class multiply_mat_vari<double, 1, CA, TB, 1> : public vari {
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
        multiply_mat_vari(const Eigen::Matrix<double, 1, CA>& A,
                          const Eigen::Matrix<TB, CA, 1>& B)
          : vari(0.0),
          size_(A.cols()),
          Ad_(ChainableStack::memalloc_.alloc_array<double>(size_)),
          Bd_(ChainableStack::memalloc_.alloc_array<double>(size_)),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {
            using Eigen::Map;
            using Eigen::VectorXd;
            using Eigen::RowVectorXd;
            for (size_type i = 0; i < size_; ++i)
              Ad_[i] = A.coeffRef(i);
            for (size_type i = 0; i < size_; ++i) {
              variRefB_[i] = B.coeffRef(i).vi_;
              Bd_[i] = B.coeffRef(i).val();
            }
            double AB
              = Eigen::Map<RowVectorXd>(Ad_, 1, size_)
              * Eigen::Map<VectorXd>(Bd_, size_, 1);
            variRefAB_ = new vari(AB, false);
          }

        virtual void chain() {
          using Eigen::RowVectorXd;
          using Eigen::VectorXd;
          using Eigen::Map;
          double adjAB;
          VectorXd adjB(size_);

          adjAB = variRefAB_->adj_;
          adjB = Map<RowVectorXd>(Ad_, 1, size_).transpose()
            * adjAB;
          for (size_type i = 0; i < size_; ++i)
            variRefB_[i]->adj_ += adjB(i);
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
     * @tparam TA Scalar type for matrix A
     * @tparam RA Rows for matrix A
     * @tparam CA Columns for matrix A, Rows for matrix B
     * @tparam CB Columns for matrix B
     */
    template <typename TA, int RA, int CA, int CB>
    class multiply_mat_vari<TA, RA, CA, double, CB> : public vari {
      public:
        int A_rows_, A_cols_, B_cols_, A_size_, B_size_;
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
        multiply_mat_vari(const Eigen::Matrix<TA, RA, CA>& A,
                          const Eigen::Matrix<double, CA, CB>& B)
          : vari(0.0),
          A_rows_(A.rows()), A_cols_(A.cols()),
          B_cols_(B.cols()), A_size_(A.size()), B_size_(B.size()),
          Ad_(ChainableStack::memalloc_.alloc_array<double>(A_size_)),
          Bd_(ChainableStack::memalloc_.alloc_array<double>(B_size_)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>(A_size_)),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>(A_rows_
                                                                  * B_cols_)) {
            using Eigen::Map;
            using Eigen::MatrixXd;
            for (size_type i = 0; i < A_size_; ++i) {
              variRefA_[i] = A.coeffRef(i).vi_;
              Ad_[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < B_size_; ++i) {
              Bd_[i] = B.coeffRef(i);
            }
            MatrixXd AB
              = Map<MatrixXd>(Ad_, A_rows_, A_cols_)
              * Map<MatrixXd>(Bd_, A_cols_, B_cols_);
            for (size_type i = 0; i < AB.size(); ++i)
              variRefAB_[i] = new vari(AB.coeffRef(i), false);
          }

        virtual void chain() {
          using Eigen::MatrixXd;
          using Eigen::Map;
          MatrixXd adjAB(A_rows_, B_cols_);
          MatrixXd adjA(A_rows_, A_cols_);

          for (size_type i = 0; i < adjAB.size(); ++i)
            adjAB(i) = variRefAB_[i]->adj_;
          adjA = adjAB * Map<MatrixXd>(Bd_, A_cols_, B_cols_).transpose();
          for (size_type i = 0; i < A_size_; ++i)
            variRefA_[i]->adj_ += adjA(i);
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
     * @tparam TA Scalar type for matrix A
     * @tparam RA Rows for matrix A
     * @tparam CA Columns for matrix A, Rows for matrix B
     * @tparam TB Scalar type for matrix B
     * @tparam CB Columns for matrix B
     */
    template <typename TA, int CA>
    class multiply_mat_vari<TA, 1, CA, double, 1> : public vari {
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
        multiply_mat_vari(const Eigen::Matrix<TA, 1, CA>& A,
                          const Eigen::Matrix<double, CA, 1>& B)
          : vari(0.0),
          size_(A.cols()),
          Ad_(ChainableStack::memalloc_.alloc_array<double>(size_)),
          Bd_(ChainableStack::memalloc_.alloc_array<double>(size_)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>(size_)) {
            using Eigen::Map;
            using Eigen::VectorXd;
            using Eigen::RowVectorXd;
            for (size_type i = 0; i < size_; ++i) {
              variRefA_[i] = A.coeffRef(i).vi_;
              Ad_[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < size_; ++i)
              Bd_[i] = B.coeffRef(i);
            double AB
              = Map<RowVectorXd>(Ad_, 1, size_)
              * Map<VectorXd>(Bd_, size_, 1);
            variRefAB_ = new vari(AB, false);
          }

        virtual void chain() {
          using Eigen::Map;
          using Eigen::VectorXd;
          using Eigen::RowVectorXd;
          double adjAB;
          RowVectorXd adjA(size_);

          adjAB = variRefAB_->adj_;
          adjA = adjAB
            * Map<VectorXd>(Bd_, size_, 1).transpose();
          for (size_type i = 0; i < size_; ++i)
            variRefA_[i]->adj_ += adjA(i);
        }
    };

    /**
     * Return the product of two scalars.
     * @tparam T1 scalar type of v
     * @tparam T2 scalar type of c
     * @param[in] v First scalar
     * @param[in] c Specified scalar
     * @return Product of scalars
     */
    template <typename T1, typename T2>
    inline typename
    boost::enable_if_c<
    (boost::is_scalar<T1>::value || boost::is_same<T1, var>::value)
    && (boost::is_scalar<T2>::value || boost::is_same<T2, var>::value),
    typename boost::math::tools::promote_args<T1, T2>::type>::type
      multiply(const T1& v, const T2& c) {
        return v * c;
      }

    /**
     * Return the product of scalar and matrix.
     * @tparam T1 scalar type v
     * @tparam T2 scalar type matrix m
     * @tparam R2 Rows matrix m
     * @tparam C2 Columns matrix m
     * @param[in] c Specified scalar
     * @param[in] m Matrix
     * @return Product of scalar and matrix
     */
    template<typename T1, typename T2, int R2, int C2>
    inline Eigen::Matrix<var, R2, C2>
    multiply(const T1& c, const Eigen::Matrix<T2, R2, C2>& m) {
      // TODO(trangucci) pull out to eliminate overpromotion of one side
      // move to matrix.hpp w. promotion?
      return to_var(m) * to_var(c);
    }

    /**
     * Return the product of scalar and matrix.
     * @tparam T1 scalar type matrix m 
     * @tparam T2 scalar type v 
     * @tparam R1 Rows matrix m
     * @tparam C1 Columns matrix m
     * @param[in] c Specified scalar
     * @param[in] m Matrix
     * @return Product of scalar and matrix
     */
    template<typename T1, int R1, int C1, typename T2>
    inline Eigen::Matrix<var, R1, C1>
    multiply(const Eigen::Matrix<T1, R1, C1>& m, const T2& c) {
      // TODO(trangucci) pull out to eliminate overpromotion of one side
      // move to matrix.hpp w. promotion?
      return to_var(m) * to_var(c);
    }

    /**
     * Return the product of two matrices.
     * @tparam TA scalar type matrix A 
     * @tparam RA Rows matrix A
     * @tparam CA Columns matrix A
     * @tparam TB scalar type matrix B 
     * @tparam RB Rows matrix B
     * @tparam CB Columns matrix B
     * @param[in] A Matrix
     * @param[in] B Matrix
     * @return Product of scalar and matrix.
     */
    template <typename TA, int RA, int CA, typename TB, int CB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value
                       || boost::is_same<TB, var>::value,
                       Eigen::Matrix<var, RA, CB> >::type
    multiply(const Eigen::Matrix<TA, RA, CA> &A,
             const Eigen::Matrix<TB, CA, CB> &B) {
      check_multiplicable("multiply", "A", A, "B", B);
      check_not_nan("multiply", "A", A);
      check_not_nan("multiply", "B", B);

      // Memory managed with the arena allocator.
      multiply_mat_vari<TA, RA, CA, TB, CB> *baseVari
        = new multiply_mat_vari<TA, RA, CA, TB, CB>(A, B);
      Eigen::Matrix<var, RA, CB> AB_v(A.rows(), B.cols());
      for (size_type i = 0; i < AB_v.size(); ++i) {
        AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
      }
      return AB_v;
    }

    /**
     * Return the scalar product of a row vector and 
     * a vector.
     * @tparam TA scalar type row vector A 
     * @tparam CA Columns matrix A
     * @tparam TB scalar type vector B 
     * @param[in] A Row vector
     * @param[in] B Column vector
     * @return Scalar product of row vector and vector
     */
    template <typename TA, int CA, typename TB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value
                       || boost::is_same<TB, var>::value, var>::type
    multiply(const Eigen::Matrix<TA, 1, CA> &A,
             const Eigen::Matrix<TB, CA, 1> &B) {
      check_multiplicable("multiply", "A", A, "B", B);
      check_not_nan("multiply", "A", A);
      check_not_nan("multiply", "B", B);

      // Memory managed with the arena allocator.
      multiply_mat_vari<TA, 1, CA, TB, 1> *baseVari
        = new multiply_mat_vari<TA, 1, CA, TB, 1>(A, B);
      var AB_v;
      AB_v.vi_ = baseVari->variRefAB_;
      return AB_v;
    }
  }
}
#endif
