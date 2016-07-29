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

namespace stan {
  namespace math {

    /**
     * Return the product of two scalars.
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
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
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
              = Eigen::Map<Eigen::Matrix<double,RA_,CA_> >(Ad_a,RA, CA) 
              * Eigen::Map<Eigen::Matrix<double,CA_,CB_> >(Bd_a,CA,CB);
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        using Eigen::Map;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, RA_, CA_> adjA(RA, CA);
        Matrix<double, CA_, CB_> adjB(CA, CB);

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Map<Matrix<double, CA_, CB_> >(Bd_a,CA,CB).transpose();
        adjB = Map<Matrix<double, RA_, CA_> >(Ad_a,RA,CA).transpose() * adjAB;
        for (size_type i = 0; i < Asize; ++i) 
            variRefA_[i]->adj_ += adjA(i); 
        for (size_type i = 0; i < Bsize; ++i) 
            variRefB_[i]->adj_ += adjB(i); 
      }
    };

    template <typename TA, int RA_, int CA_, typename TB>
    class multiply_mat_scal_vari : public vari {
    public:
      int RA, CA, Asize;  // A.cols() = B.rows()
      const double Bd_;
      double* Ad_a;
      vari** variRefA_;
      vari* variRefB_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_scal_vari(const Eigen::Matrix<TA, RA_, CA_>& A,
                            const TB& B)
        : vari(0.0),
          RA(A.rows()), CA(A.cols()), Asize(A.size()),
          Bd_(value_of(B)),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())) {
            for (size_type i = 0; i < A.size(); ++i) { 
                variRefA_[i] = A.coeffRef(i).vi_;
                Ad_a[i] = A.coeffRef(i).val();
            }
            Eigen::Matrix<double, RA_, CA_> AB
              = Eigen::Map<Eigen::Matrix<double,RA_,CA_> >(Ad_a,RA,CA) * Bd_;
            variRefB_ = B.vi_;
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Map;
        using Eigen::Infinity;
        Matrix<double, RA_, CA_> adjAB(RA, CA);
        Matrix<double, RA_, CA_> adjA(RA, CA);
        double adjB;

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Bd_;
        adjB = (Map<Eigen::Array<double,RA_,CA_> >(Ad_a,RA,CA) * adjAB.array()).sum();
        for (size_type i = 0; i < Asize; ++i) 
            variRefA_[i]->adj_ += adjA(i); 
        variRefB_->adj_ += adjB; 
      }
    };

    template <int RA_, int CA_, typename TB>
    class multiply_mat_scal_vari<double,RA_,CA_,TB> : public vari {
    public:
      int RA, CA, Asize;  // A.cols() = B.rows()
      const double Bd_;
      double* Ad_a;
      vari* variRefB_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_scal_vari(const Eigen::Matrix<double, RA_, CA_>& A,
                            const TB& B)
        : vari(0.0),
          RA(A.rows()), CA(A.cols()), Asize(A.size()),
          Bd_(value_of(B)),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())) {
            for (size_type i = 0; i < A.size(); ++i) { 
                Ad_a[i] = A.coeffRef(i);
            }
            Eigen::Matrix<double, RA_, CA_> AB
              = Eigen::Map<Eigen::Matrix<double,RA_,CA_> >(Ad_a,RA,CA) * Bd_;
            variRefB_ = B.vi_;
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Map;
        using Eigen::Infinity;
        Matrix<double, RA_, CA_> adjAB(RA, CA);
        double adjB;

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjB = (Map<Eigen::Array<double,RA_,CA_> >(Ad_a,RA,CA) * adjAB.array()).sum();
        variRefB_->adj_ += adjB; 
      }
    };

    template <typename TA, int RA_, int CA_>
    class multiply_mat_scal_vari<TA, RA_, CA_, double> : public vari {
    public:
      int RA, CA, Asize;  // A.cols() = B.rows()
      const double Bd_;
      double* Ad_a;
      vari** variRefA_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
       * @param matrix A
       * @param matrix B
       * */
     multiply_mat_scal_vari(const Eigen::Matrix<TA, RA_, CA_>& A,
                            const double B)
        : vari(0.0),
          RA(A.rows()), CA(A.cols()), Asize(A.size()),
          Bd_(value_of(B)),
          Ad_a(ChainableStack::memalloc_.alloc_array<double>
               (A.rows() * A.cols())),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())) {
            for (size_type i = 0; i < A.size(); ++i) { 
                variRefA_[i] = A.coeffRef(i).vi_;
                Ad_a[i] = A.coeffRef(i).val();
            }
            Eigen::Matrix<double, RA_, CA_> AB
              = Eigen::Map<Eigen::Matrix<double,RA_,CA_> >(Ad_a,RA,CA) * Bd_;
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        Matrix<double, RA_, CA_> adjAB(RA, CA);
        Matrix<double, RA_, CA_> adjA(RA, CA);

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Bd_;
        for (size_type i = 0; i < Asize; ++i) 
            variRefA_[i]->adj_ += adjA(i); 
      }
    };

    template <typename TA, int CA_, typename TB>
    class multiply_mat_vari<TA, 1, CA_, TB, 1> : public vari {
    public:
      int CA;  // A.cols() = B.rows()
      double* Ad_a; // row_vector
      double* Bd_a; // column_vector
      vari** variRefA_;
      vari** variRefB_;
      vari* variRefAB_;

      /* ctor for multiply function
       *
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
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
            double AB = Eigen::Map<Eigen::Matrix<double, 1, CA_> > (Ad_a,1,CA)
                        * Eigen::Map<Eigen::Matrix<double, CA_, 1> > (Bd_a,CA,1);
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Map;
        using Eigen::Infinity;
        double adjAB;
        Matrix<double, 1, CA_> adjA(1, CA);
        Matrix<double, CA_, 1> adjB(CA, 1);

        adjAB = variRefAB_->adj_;
        adjA = adjAB * Map<Matrix<double,CA_,1> >(Bd_a,CA,1).transpose();
        adjB = Map<Matrix<double,1,CA_> >(Ad_a,1,CA).transpose() * adjAB;
        for (size_type i = 0; i < CA; ++i) 
            variRefA_[i]->adj_ += adjA(i); 
        for (size_type i = 0; i < CA; ++i) 
            variRefB_[i]->adj_ += adjB(i); 
      }
    };

    template <int RA_, int CA_, typename TB, int CB_>
    class multiply_mat_vari<double,RA_,CA_,TB,CB_> : public vari {
    public:
      int RA, CA, CB, Asize, Bsize;  // A.cols() = B.rows()
      double* Ad_a;
      double* Bd_a;
      vari** variRefB_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
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
              = Eigen::Map<Eigen::Matrix<double,RA_,CA_> >(Ad_a,RA, CA) 
              * Eigen::Map<Eigen::Matrix<double,CA_,CB_> >(Bd_a,CA,CB);
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        using Eigen::Map;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, CA_, CB_> adjB(CA, CB);

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjB = Map<Matrix<double, RA_, CA_> >(Ad_a,RA,CA).transpose() * adjAB;
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
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
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
              = Eigen::Map<Eigen::Matrix<double,1,CA_> >(Ad_a,1, CA) 
              * Eigen::Map<Eigen::Matrix<double,CA_,1> >(Bd_a,CA,1);
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        double adjAB;
        Matrix<double, CA_, 1> adjB(CA, 1);

        adjAB = variRefAB_->adj_;
        adjB = Eigen::Map<Eigen::Matrix<double,1,CA_> >(Ad_a,1,CA).transpose() * adjAB;
        for (size_type i = 0; i < CA; ++i) 
            variRefB_[i]->adj_ += adjB(i); 
      }
    };

    template <typename TA, int RA_, int CA_, int CB_>
    class multiply_mat_vari<TA,RA_,CA_,double,CB_> : public vari {
    public:
      int RA, CA, CB, Asize, Bsize;  // A.cols() = B.rows()
      double* Ad_a;
      double* Bd_a;
      vari** variRefA_;
      vari** variRefAB_;

      /* ctor for multiply function
       *
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
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
              = Eigen::Map<Eigen::Matrix<double,RA_,CA_> >(Ad_a,RA, CA) 
              * Eigen::Map<Eigen::Matrix<double,CA_,CB_> >(Bd_a,CA,CB);
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        using Eigen::Map;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, RA_, CA_> adjA(RA, CA);

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Map<Matrix<double, CA_, CB_> >(Bd_a,CA,CB).transpose();
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
       * Stores varis for A
       * Stores varis for B
       * Stores varis for AB
       *
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
              = Eigen::Map<Eigen::Matrix<double,1,CA_> >(Ad_a,1, CA) 
              * Eigen::Map<Eigen::Matrix<double,CA_,1> >(Bd_a,CA,1);
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        double adjAB;
        Matrix<double, 1, CA_> adjA(1, CA);

        adjAB = variRefAB_->adj_;
        adjA = adjAB * Eigen::Map<Eigen::Matrix<double,CA_,1> >(Bd_a,CA,1).transpose();
        for (size_type i = 0; i < CA; ++i) 
            variRefA_[i]->adj_ += adjA(i); 
      }
    };

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
      stan::math::check_not_nan("multiply","A", A);
      stan::math::check_not_nan("multiply","B", B);

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      multiply_mat_vari<TA,RA,CA,TB,CB> *baseVari
        = new multiply_mat_vari<TA,RA,CA,TB,CB>(A, B);
      Eigen::Matrix<var, RA, CB> AB_v(A.rows(), B.cols());
      for (size_type i = 0; i < AB_v.size(); ++i) {
          AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
        }
      return AB_v;
    }

    template <typename TA, int RA, int CA, typename TB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value ||
                        boost::is_same<TB, var>::value,
                        Eigen::Matrix<var, RA, CA> >::type
    multiply(const Eigen::Matrix<TA, RA, CA>& A,
             const TB& B) {
      stan::math::check_not_nan("multiply","A", A);
      stan::math::check_not_nan("multiply","B", B);

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      multiply_mat_scal_vari<TA,RA,CA,TB> *baseVari
        = new multiply_mat_scal_vari<TA,RA,CA,TB>(A, B);
      Eigen::Matrix<var, RA, CA> AB_v(A.rows(), A.cols());
      for (size_type i = 0; i < AB_v.size(); ++i) {
          AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
        }
      return AB_v;
    }

    template <typename TA, int RA, int CA, typename TB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value ||
                        boost::is_same<TB, var>::value,
                        Eigen::Matrix<var, RA, CA> >::type
    multiply(const TB& B,
             const Eigen::Matrix<TA, RA, CA>& A) {
      stan::math::check_not_nan("multiply","A", A);
      stan::math::check_not_nan("multiply","B", B);

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      multiply_mat_scal_vari<TA,RA,CA,TB> *baseVari
        = new multiply_mat_scal_vari<TA,RA,CA,TB>(A, B);
      Eigen::Matrix<var, RA, CA> AB_v(A.rows(), A.cols());
      for (size_type i = 0; i < AB_v.size(); ++i) {
          AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
        }
      return AB_v;
    }

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
      stan::math::check_not_nan("multiply","A", A);
      stan::math::check_not_nan("multiply","B", B);

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      multiply_mat_vari<TA,1,CA,TB,1> *baseVari
        = new multiply_mat_vari<TA,1,CA,TB,1>(A, B);
      var AB_v;
      AB_v.vi_ = baseVari->variRefAB_;
      return AB_v;
    }
  }
}
#endif
