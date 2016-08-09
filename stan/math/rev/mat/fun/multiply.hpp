#ifndef STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>

namespace stan {
  namespace math {

    template <typename TA, int RA_, int CA_, typename TB, int CB_>
    class multiply_mat_vari : public vari {
    public:
      int RA, CA, CB, Asize, Bsize;  // A.cols() = B.rows()
      double* Ad_a;
      double* Bd_a;
      vari** variRefA_;
      vari** variRefB_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores pointers to varis for A
       * Stores pointers to varis for B
       * Stores pointers to varis for AB
       * All memory allocated on ChainableStack
       *
       * @tparam TA Scalar type for matrix A
       * @tparam RA_ Rows for matrix A
       * @tparam CA_ Columns for matrix A, Rows for matrix B
       * @tparam TB Scalar type for matrix B
       * @tparam CB_ Columns for matrix B
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_vari(const Eigen::Matrix<TA, RA_, CA_>& A,
                       const Eigen::Matrix<TB, CA_, CB_>& B)
        : vari(0.0),
          RA(A.rows()), CA(A.cols()),
          CB(B.cols()), Asize(A.size()), Bsize(B.size()),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          Bd_a(ChainableStack::memalloc_.alloc_array<double>
               (B.rows() * B.cols())),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * B.cols())) {
            for (size_type i = 0; i < A.size(); ++i) {
                variRefA_[i] = A.coeffRef(i).vi_;
                Ad_a[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < B.size(); ++i) {
                variRefB_[i] = B.coeffRef(i).vi_;
                Bd_a[i] = B.coeffRef(i).val();
            }
            Eigen::Matrix<double, RA_, CB_> AB
              = Eigen::Map<Eigen::Matrix<double, RA_, CA_> >(Ad_a, RA, CA)
              * Eigen::Map<Eigen::Matrix<double, CA_, CB_> >(Bd_a, CA, CB);
            for (size_type i = 0; i < AB.size(); ++i)
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Map;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, RA_, CA_> adjA(RA, CA);
        Matrix<double, CA_, CB_> adjB(CA, CB);

        for (size_type i = 0; i < adjAB.size(); ++i)
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Map<Matrix<double, CA_, CB_> >(Bd_a, CA, CB).transpose();
        adjB = Map<Matrix<double, RA_, CA_> >(Ad_a, RA, CA).transpose() * adjAB;
        for (size_type i = 0; i < Asize; ++i)
            variRefA_[i]->adj_ += adjA(i);
        for (size_type i = 0; i < Bsize; ++i)
            variRefB_[i]->adj_ += adjB(i);
      }
    };

    template <typename TA, int CA_, typename TB>
    class multiply_mat_vari<TA, 1, CA_, TB, 1> : public vari {
    public:
      int CA;
      double* Ad_a;
      double* Bd_a;
      vari** variRefA_;
      vari** variRefB_;
      vari* variRefAB_;

      /* ctor for multiply function
       *
       * Stores pointers to varis for A
       * Stores pointers to varis for B
       * Stores pointers to varis for AB
       * All memory allocated on ChainableStack
       *
       * @tparam TA Scalar type for matrix A
       * @tparam CA_ Columns for matrix A, Rows for matrix B
       * @tparam TB Scalar type for matrix B
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_vari(const Eigen::Matrix<TA, 1, CA_>& A,
                       const Eigen::Matrix<TB, CA_, 1>& B)
        : vari(0.0),
          CA(A.cols()),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.cols())),
          Bd_a(ChainableStack::memalloc_.alloc_array<double>
               (B.rows())),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())) {
            for (size_type i = 0; i < A.size(); ++i) {
                variRefA_[i] = A.coeffRef(i).vi_;
                Ad_a[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < B.size(); ++i) {
                variRefB_[i] = B.coeffRef(i).vi_;
                Bd_a[i] = B.coeffRef(i).val();
            }
            double AB = Eigen::Map<
              Eigen::Matrix<double, 1, CA_> > (Ad_a, 1, CA)
                        * Eigen::Map<
              Eigen::Matrix<double, CA_, 1> > (Bd_a, CA, 1);
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Map;
        double adjAB;
        Matrix<double, 1, CA_> adjA(1, CA);
        Matrix<double, CA_, 1> adjB(CA, 1);

        adjAB = variRefAB_->adj_;
        adjA = adjAB * Map<Matrix<double, CA_, 1> >(Bd_a, CA, 1).transpose();
        adjB = Map<Matrix<double, 1, CA_> >(Ad_a, 1, CA).transpose() * adjAB;
        for (size_type i = 0; i < CA; ++i)
            variRefA_[i]->adj_ += adjA(i);
        for (size_type i = 0; i < CA; ++i)
            variRefB_[i]->adj_ += adjB(i);
      }
    };

    template <int RA_, int CA_, typename TB, int CB_>
    class multiply_mat_vari<double, RA_, CA_, TB, CB_> : public vari {
    public:
      int RA, CA, CB, Asize, Bsize;  // A.cols() = B.rows()
      double* Ad_a;
      double* Bd_a;
      vari** variRefB_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores pointers to varis for B
       * Stores pointers to varis for AB
       * All memory allocated on ChainableStack
       *
       * @tparam RA_ Rows for matrix A
       * @tparam CA_ Columns for matrix A, Rows for matrix B
       * @tparam TB Scalar type for matrix B
       * @tparam CB_ Columns for matrix B
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_vari(const Eigen::Matrix<double, RA_, CA_>& A,
                       const Eigen::Matrix<TB, CA_, CB_>& B)
        : vari(0.0),
          RA(A.rows()), CA(A.cols()),
          CB(B.cols()), Asize(A.size()), Bsize(B.size()),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          Bd_a(ChainableStack::memalloc_.alloc_array<double>
               (B.rows() * B.cols())),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * B.cols())) {
            for (size_type i = 0; i < A.size(); ++i) {
                Ad_a[i] = A.coeffRef(i);
            }
            for (size_type i = 0; i < B.size(); ++i) {
                variRefB_[i] = B.coeffRef(i).vi_;
                Bd_a[i] = B.coeffRef(i).val();
            }
            Eigen::Matrix<double, RA_, CB_> AB
              = Eigen::Map<Eigen::Matrix<double, RA_, CA_> >(Ad_a, RA, CA)
              * Eigen::Map<Eigen::Matrix<double, CA_, CB_> >(Bd_a, CA, CB);
            for (size_type i = 0; i < AB.size(); ++i)
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Map;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, CA_, CB_> adjB(CA, CB);

        for (size_type i = 0; i < adjAB.size(); ++i)
            adjAB(i) = variRefAB_[i]->adj_;
        adjB = Map<Matrix<double, RA_, CA_> >(Ad_a, RA, CA).transpose() * adjAB;
        for (size_type i = 0; i < Bsize; ++i)
            variRefB_[i]->adj_ += adjB(i);
      }
    };

    template <int CA_, typename TB>
    class multiply_mat_vari<double, 1, CA_, TB, 1> : public vari {
    public:
      int CA;  // A.cols() = B.rows()
      double* Ad_a;
      double* Bd_a;
      vari** variRefB_;
      vari* variRefAB_;

      /* ctor for multiply function
       *
       * Stores pointers to varis for B
       * Stores pointers to varis for AB
       * All memory allocated on ChainableStack
       *
       * @tparam CA_ Columns for matrix A, Rows for matrix B
       * @tparam TB Scalar type for matrix B
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_vari(const Eigen::Matrix<double, 1, CA_>& A,
                       const Eigen::Matrix<TB, CA_, 1>& B)
        : vari(0.0),
          CA(A.cols()),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          Bd_a(ChainableStack::memalloc_.alloc_array<double>
               (B.rows() * B.cols())),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())) {
            for (size_type i = 0; i < A.size(); ++i) {
                Ad_a[i] = A.coeffRef(i);
            }
            for (size_type i = 0; i < B.size(); ++i) {
                variRefB_[i] = B.coeffRef(i).vi_;
                Bd_a[i] = B.coeffRef(i).val();
            }
            double AB
              = Eigen::Map<Eigen::Matrix<double, 1, CA_> >(Ad_a, 1, CA)
              * Eigen::Map<Eigen::Matrix<double, CA_, 1> >(Bd_a, CA, 1);
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        double adjAB;
        Matrix<double, CA_, 1> adjB(CA, 1);

        adjAB = variRefAB_->adj_;
        adjB = Eigen::Map<
          Eigen::Matrix<double, 1, CA_> >(Ad_a, 1, CA).transpose()
          * adjAB;
        for (size_type i = 0; i < CA; ++i)
            variRefB_[i]->adj_ += adjB(i);
      }
    };

    template <typename TA, int RA_, int CA_, int CB_>
    class multiply_mat_vari<TA, RA_, CA_, double, CB_> : public vari {
    public:
      int RA, CA, CB, Asize, Bsize;  // A.cols() = B.rows()
      double* Ad_a;
      double* Bd_a;
      vari** variRefA_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores pointers to varis for A
       * Stores pointers to varis for AB
       * All memory allocated on ChainableStack
       *
       * @tparam TA Scalar type for matrix A
       * @tparam RA_ Rows for matrix A
       * @tparam CA_ Columns for matrix A, Rows for matrix B
       * @tparam CB_ Columns for matrix B
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_vari(const Eigen::Matrix<TA, RA_, CA_>& A,
                       const Eigen::Matrix<double, CA_, CB_>& B)
        : vari(0.0),
          RA(A.rows()), CA(A.cols()),
          CB(B.cols()), Asize(A.size()), Bsize(B.size()),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          Bd_a(ChainableStack::memalloc_.alloc_array<double>
               (B.rows() * B.cols())),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * B.cols())) {
            for (size_type i = 0; i < A.size(); ++i) {
                variRefA_[i] = A.coeffRef(i).vi_;
                Ad_a[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < B.size(); ++i) {
                Bd_a[i] = B.coeffRef(i);
            }
            Eigen::Matrix<double, RA_, CB_> AB
              = Eigen::Map<Eigen::Matrix<double, RA_, CA_> >(Ad_a, RA, CA)
              * Eigen::Map<Eigen::Matrix<double, CA_, CB_> >(Bd_a, CA, CB);
            for (size_type i = 0; i < AB.size(); ++i)
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Map;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, RA_, CA_> adjA(RA, CA);

        for (size_type i = 0; i < adjAB.size(); ++i)
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Map<Matrix<double, CA_, CB_> >(Bd_a, CA, CB).transpose();
        for (size_type i = 0; i < Asize; ++i)
            variRefA_[i]->adj_ += adjA(i);
      }
    };

    template <typename TA, int CA_>
    class multiply_mat_vari<TA, 1, CA_, double, 1> : public vari {
    public:
      int CA;  // A.cols() = B.rows()
      double* Ad_a;
      double* Bd_a;
      vari** variRefA_;
      vari* variRefAB_;

      /* ctor for multiply function
       *
       * Stores pointers to varis for A
       * Stores pointers to varis for AB
       * All memory allocated on ChainableStack
       *
       * @tparam TA Scalar type for matrix A
       * @tparam CA_ Columns for matrix A, Rows for matrix B
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_vari(const Eigen::Matrix<TA, 1, CA_>& A,
                       const Eigen::Matrix<double, CA_, 1>& B)
        : vari(0.0),
          CA(A.cols()),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          Bd_a(ChainableStack::memalloc_.alloc_array<double>
               (B.rows() * B.cols())),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())) {
            for (size_type i = 0; i < A.size(); ++i) {
                variRefA_[i] = A.coeffRef(i).vi_;
                Ad_a[i] = A.coeffRef(i).val();
            }
            for (size_type i = 0; i < B.size(); ++i) {
                Bd_a[i] = B.coeffRef(i);
            }
            double AB
              = Eigen::Map<Eigen::Matrix<double, 1, CA_> >(Ad_a, 1, CA)
              * Eigen::Map<Eigen::Matrix<double, CA_, 1> >(Bd_a, CA, 1);
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        double adjAB;
        Matrix<double, 1, CA_> adjA(1, CA);

        adjAB = variRefAB_->adj_;
        adjA = adjAB
          * Eigen::Map<Eigen::Matrix<double, CA_, 1> >(Bd_a, CA, 1).transpose();
        for (size_type i = 0; i < CA; ++i)
            variRefA_[i]->adj_ += adjA(i);
      }
    };

    /**
     * Return the product of two scalars.
     * @tparam T1 scalar type of v
     * @tparam T2 scalar type of c
     * @param[in] v First scalar.
     * @param[in] c Specified scalar.
     * @return Product of scalars.
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
     * @param[in] c Specified scalar.
     * @param[in] m Matrix.
     * @return Product of scalar and matrix.
     */
    template<typename T1, typename T2, int R2, int C2>
    inline Eigen::Matrix<var, R2, C2>
    multiply(const T1& c, const Eigen::Matrix<T2, R2, C2>& m) {
      // FIXME:  pull out to eliminate overpromotion of one side
      // move to matrix.hpp w. promotion?
      return to_var(m) * to_var(c);
    }

    /**
     * Return the product of scalar and matrix.
     * @tparam T1 scalar type matrix m 
     * @tparam T2 scalar type v 
     * @tparam R1 Rows matrix m
     * @tparam C1 Columns matrix m
     * @param[in] c Specified scalar.
     * @param[in] m Matrix.
     * @return Product of scalar and matrix.
     */
    template<typename T1, int R1, int C1, typename T2>
    inline Eigen::Matrix<var, R1, C1>
    multiply(const Eigen::Matrix<T1, R1, C1>& m, const T2& c) {
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
     * @param[in] A Matrix.
     * @param[in] B Matrix.
     * @return Product of scalar and matrix.
     */
    template <typename TA, int RA, int CA, typename TB, int CB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value ||
                        boost::is_same<TB, var>::value,
                        Eigen::Matrix<var, RA, CB> >::type
    multiply(const Eigen::Matrix<TA, RA, CA> &A,
             const Eigen::Matrix<TB, CA, CB> &B) {
      stan::math::check_multiplicable("multiply",
                                                "A", A,
                                                "B", B);
      stan::math::check_not_nan("multiply", "A", A);
      stan::math::check_not_nan("multiply", "B", B);

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      multiply_mat_vari<TA, RA, CA, TB, CB> *baseVari
        = new multiply_mat_vari<TA, RA, CA, TB, CB>(A, B);
      Eigen::Matrix<var, RA, CB> AB_v(A.rows(), B.cols());
      for (size_type i = 0; i < AB_v.size(); ++i) {
          AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
        }
      return AB_v;
    }

    /**
     * Return the scalar product of two vectors.
     * @tparam TA scalar type row vector A 
     * @tparam CA Columns matrix A
     * @tparam TB scalar type vector B 
     * @param[in] A Row vector.
     * @param[in] B Column vector.
     * @return Scalar product of row vector and vector.
     */
    template <typename TA, int CA, typename TB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value ||
                       boost::is_same<TB, var>::value,
                       var>::type
    multiply(const Eigen::Matrix<TA, 1, CA> &A,
             const Eigen::Matrix<TB, CA, 1> &B) {
      stan::math::check_multiplicable("multiply",
                                                "A", A,
                                                "B", B);
      stan::math::check_not_nan("multiply", "A", A);
      stan::math::check_not_nan("multiply", "B", B);

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      multiply_mat_vari<TA, 1, CA, TB, 1> *baseVari
        = new multiply_mat_vari<TA, 1, CA, TB, 1>(A, B);
      var AB_v;
      AB_v.vi_ = baseVari->variRefAB_;
      return AB_v;
    }
  }
}
#endif
