#ifndef STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_MULTIPLY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>

namespace stan {
  namespace math {

    class multiply_mat_vv_vari : public vari {
    public:
      int RA, CA, RB, CB;  // A.rows() = A.cols()
      const Eigen::Matrix<double, -1, -1> Ad_;
      const Eigen::Matrix<double, -1, -1> Bd_;
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
     multiply_mat_vv_vari(const Eigen::Matrix<var, -1, -1>& A,
                          const Eigen::Matrix<var, -1, -1>& B)
        : vari(0.0),
          RA(A.rows()), CA(A.cols()), 
          RB(B.rows()), CB(B.cols()),
          Ad_(value_of(A)), Bd_(value_of(B)),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * A.cols())),
          variRefB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (B.rows() * B.cols())),
          variRefAB_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * B.cols())) {
            Eigen::Matrix<double, -1, -1> AB
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
        Matrix<double, -1, -1> adjAB(RA, CB);
        Matrix<double, -1, -1> adjA(RA, CA);
        Matrix<double, -1, -1> adjB(RB, CB);

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

    Eigen::Matrix<var, -1, -1>
    multiply(const Eigen::Matrix<var, -1, -1> &A,
             const Eigen::Matrix<var, -1, -1> &B) {
      stan::math::check_multiplicable("multiply",
                                                "A", A,
                                                "B", B);
      stan::math::check_not_nan("multiply","A", A);
      stan::math::check_not_nan("multiply","B", B);

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      multiply_mat_vv_vari *baseVari
        = new multiply_mat_vv_vari(A, B);
      Eigen::Matrix<var, -1, -1> AB_v(A.rows(), B.cols());
      for (size_type i = 0; i < AB_v.size(); ++i) {
          AB_v.coeffRef(i).vi_ = baseVari->variRefAB_[i];
        }
      return AB_v;
    }
  }
}
#endif
