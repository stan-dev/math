#ifndef STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>

namespace stan {
  namespace math {

    template <typename TA, int RA_, int CA_, typename TB, int CB_>
    class multiply_mat_vari : public vari {
    public:
      int RA, CA, CB;  // A.cols() = B.rows()
      const Eigen::Matrix<double, RA_, CA_> Ad_;
      const Eigen::Matrix<double, CA_, CB_> Bd_;
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
          CB(B.cols()),
          Ad_(value_of(A)), Bd_(value_of(B)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * B.cols())) {
            Eigen::Matrix<double, RA_, CB_> AB
              = Ad_ * Bd_;
            for (size_type i = 0; i < A.size(); ++i) 
                variRefA_[i] = A.coeffRef(i).vi_;
            for (size_type i = 0; i < B.size(); ++i) 
                variRefB_[i] = B.coeffRef(i).vi_;
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, RA_, CA_> adjA(RA, CA);
        Matrix<double, CA_, CB_> adjB(CA, CB);

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Bd_.transpose();
        adjB = Ad_.transpose() * adjAB;
        for (size_type i = 0; i < Ad_.size(); ++i) 
            variRefA_[i]->adj_ += adjA(i); 
        for (size_type i = 0; i < Bd_.size(); ++i) 
            variRefB_[i]->adj_ += adjB(i); 
      }
    };

    template <typename TA, int CA_, typename TB>
    class multiply_mat_vari<TA, 1, CA_, TB, 1> : public vari {
    public:
      int RA, CA, CB;  // A.cols() = B.rows()
      const Eigen::Matrix<double, 1, CA_> Ad_;
      const Eigen::Matrix<double, CA_, 1> Bd_;
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
          RA(A.rows()), CA(A.cols()), 
          CB(B.cols()),
          Ad_(value_of(A)), Bd_(value_of(B)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())) {
            double AB = Ad_ * Bd_;
            for (size_type i = 0; i < A.size(); ++i) 
                variRefA_[i] = A.coeffRef(i).vi_;
            for (size_type i = 0; i < B.size(); ++i) 
                variRefB_[i] = B.coeffRef(i).vi_;
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        double adjAB;
        Matrix<double, 1, CA_> adjA(RA, CA);
        Matrix<double, CA_, 1> adjB(CA, CB);

        adjAB = variRefAB_->adj_;
        adjA = adjAB * Bd_.transpose();
        adjB = Ad_.transpose() * adjAB;
        for (size_type i = 0; i < Ad_.size(); ++i) 
            variRefA_[i]->adj_ += adjA(i); 
        for (size_type i = 0; i < Bd_.size(); ++i) 
            variRefB_[i]->adj_ += adjB(i); 
      }
    };

    template <int RA_, int CA_, typename TB, int CB_>
    class multiply_mat_vari<double,RA_,CA_,TB,CB_> : public vari {
    public:
      int RA, CA, CB;  // A.cols() = B.rows()
      const Eigen::Matrix<double, RA_, CA_> Ad_;
      const Eigen::Matrix<double, CA_, CB_> Bd_;
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
          CB(B.cols()),
          Ad_(A), Bd_(value_of(B)),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * B.cols())) {
            Eigen::Matrix<double, RA_, CB_> AB
              = Ad_ * Bd_;
            for (size_type i = 0; i < B.size(); ++i) 
                variRefB_[i] = B.coeffRef(i).vi_;
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, CA_, CB_> adjB(CA, CB);

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjB = Ad_.transpose() * adjAB;
        for (size_type i = 0; i < Bd_.size(); ++i) 
            variRefB_[i]->adj_ += adjB(i); 
      }
    };

    template <int CA_, typename TB>
    class multiply_mat_vari<double, 1, CA_, TB, 1> : public vari {
    public:
      int RA, CA, CB;  // A.cols() = B.rows()
      const Eigen::Matrix<double, 1, CA_> Ad_;
      const Eigen::Matrix<double, CA_, 1> Bd_;
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
          RA(A.rows()), CA(A.cols()), 
          CB(B.cols()),
          Ad_(value_of(A)), Bd_(value_of(B)),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())) {
            double AB = Ad_ * Bd_;
            for (size_type i = 0; i < B.size(); ++i) 
                variRefB_[i] = B.coeffRef(i).vi_;
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        double adjAB;
        Matrix<double, CA_, 1> adjB(CA, CB);

        adjAB = variRefAB_->adj_;
        adjB = Ad_.transpose() * adjAB;
        for (size_type i = 0; i < Bd_.size(); ++i) 
            variRefB_[i]->adj_ += adjB(i); 
      }
    };

    template <typename TA, int RA_, int CA_, int CB_>
    class multiply_mat_vari<TA,RA_,CA_,double,CB_> : public vari {
    public:
      int RA, CA, CB;  // A.cols() = B.rows()
      const Eigen::Matrix<double, RA_, CA_> Ad_;
      const Eigen::Matrix<double, CA_, CB_> Bd_;
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
          CB(B.cols()),
          Ad_(value_of(A)), Bd_(value_of(B)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * B.cols())) {
            Eigen::Matrix<double, RA_, CB_> AB
              = Ad_ * Bd_;
            for (size_type i = 0; i < A.size(); ++i) 
                variRefA_[i] = A.coeffRef(i).vi_;
            for (size_type i = 0; i < AB.size(); ++i) 
                variRefAB_[i] = new vari(AB.coeffRef(i), false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        Matrix<double, RA_, CB_> adjAB(RA, CB);
        Matrix<double, RA_, CA_> adjA(RA, CA);

        for (size_type i = 0; i < adjAB.size(); ++i) 
            adjAB(i) = variRefAB_[i]->adj_;
        adjA = adjAB * Bd_.transpose();
        for (size_type i = 0; i < Ad_.size(); ++i) 
            variRefA_[i]->adj_ += adjA(i); 
      }
    };

    template <typename TA, int CA_>
    class multiply_mat_vari<TA, 1, CA_, double, 1> : public vari {
    public:
      int RA, CA, CB;  // A.cols() = B.rows()
      const Eigen::Matrix<double, 1, CA_> Ad_;
      const Eigen::Matrix<double, CA_, 1> Bd_;
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
          RA(A.rows()), CA(A.cols()), 
          CB(B.cols()),
          Ad_(value_of(A)), Bd_(value_of(B)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())) {
            double AB = Ad_ * Bd_;
            for (size_type i = 0; i < A.size(); ++i) 
                variRefA_[i] = A.coeffRef(i).vi_;
            variRefAB_ = new vari(AB, false);
        }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Infinity;
        double adjAB;
        Matrix<double, 1, CA_> adjA(RA, CA);

        adjAB = variRefAB_->adj_;
        adjA = adjAB * Bd_.transpose();
        for (size_type i = 0; i < Ad_.size(); ++i) 
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
