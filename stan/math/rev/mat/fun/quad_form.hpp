#ifndef STAN_MATH_REV_MAT_FUN_QUAD_FORM_HPP
#define STAN_MATH_REV_MAT_FUN_QUAD_FORM_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/quad_form.hpp>

namespace stan {
  namespace math {

    template <typename T, int R, int C>
    inline Eigen::Map<Eigen::Matrix<T, R, C> >
    memalloc_matrix_map(const Eigen::Matrix<T, R, C>& x) {
      T* x_mem = ChainableStack::memalloc_.template alloc_array<T>(x.size());
      for (int i = 0; i < x.size(); ++i)
        x_mem[i] = x(i);
      return Eigen::Map<Eigen::Matrix<T, R, C> >(x_mem, x.rows(), x.cols());
    }

    template <typename TA, int RA, int CA, typename TB, int RB, int CB>
    class quad_form_vari : public vari {
    protected:
      inline void chainA(const Eigen::Matrix<double, RA, CA>& A,
                         const Eigen::Matrix<double, RB, CB>& Bd,
                         const Eigen::Matrix<double, CB, CB>& adjC) {}

      inline void chainB(const Eigen::Matrix<double, RB, CB>& B,
                         const Eigen::Matrix<double, RA, CA>& Ad,
                         const Eigen::Matrix<double, RB, CB>& Bd,
                         const Eigen::Matrix<double, CB, CB>& adjC) {}

      inline void chainA(const Eigen::Matrix<var, RA, CA>& A,
                         const Eigen::Matrix<double, RB, CB>& Bd,
                         const Eigen::Matrix<double, CB, CB>& adjC) {
        Eigen::Matrix<double, RA, CA> adjA(Bd * adjC * Bd.transpose());
        for (int j = 0; j < A.cols(); j++) {
          for (int i = 0; i < A.rows(); i++) {
            A(i, j).vi_->adj_ += adjA(i, j);
          }
        }
      }

      inline void chainB(const Eigen::Matrix<var, RB, CB>& B,
                         const Eigen::Matrix<double, RA, CA>& Ad,
                         const Eigen::Matrix<double, RB, CB>& Bd,
                         const Eigen::Matrix<double, CB, CB>& adjC) {
        Eigen::Matrix<double, RA, CA> adjB(Ad * Bd * adjC.transpose()
                                           + Ad.transpose()*Bd*adjC);
        for (int j = 0; j < B.cols(); j++)
          for (int i = 0; i < B.rows(); i++)
            B(i, j).vi_->adj_ += adjB(i, j);
      }

      inline void chainAB(const Eigen::Matrix<TA, RA, CA>& A,
                          const Eigen::Matrix<TB, RB, CB>& B,
                          const Eigen::Matrix<double, RA, CA>& Ad,
                          const Eigen::Matrix<double, RB, CB>& Bd,
                          const Eigen::Matrix<double, CB, CB>& adjC) {
        chainA(A, Bd, adjC);
        chainB(B, Ad, Bd, adjC);
      }

      inline void compute(const Eigen::Matrix<double, RA, CA>& A,
                          const Eigen::Matrix<double, RB, CB>& B) {
        Eigen::Matrix<double, CB, CB> Cd(B.transpose() * A * B);
        if (sym_) {
          for (int j = 0; j < C_.cols(); j++)
            for (int i = 0; i < C_.rows(); i++)
              C_(i, j) = var(new vari(0.5 * (Cd(i, j) + Cd(j, i)), false));
        } else {
          for (int j = 0; j < C_.cols(); j++)
            for (int i = 0; i < C_.rows(); i++)
              C_(i, j) = var(new vari(Cd(i, j), false));
        }
      }

    public:
      quad_form_vari(const Eigen::Matrix<TA, RA, CA>& A,
                     const Eigen::Matrix<TB, RB, CB>& B,
                     bool sym = false)
        : vari(0.0),
          A_(A), B_(B), C_(B.cols(), B.cols()),
          map_A_(memalloc_matrix_map(A)),
          map_B_(memalloc_matrix_map(B)),
          map_C_(ChainableStack::memalloc_
                 .template alloc_array<var>(B.rows() * B.rows()),
                 B.cols(), B.cols()) {
        compute(value_of(A), value_of(B));
        for (int n = 0; n < C_.size(); ++n)
          map_C_(n) = C_(n);
      }

      virtual void chain() {
        Eigen::Matrix<double, CB, CB> adjC(C_.rows(), C_.cols());
        for (int j = 0; j < C_.cols(); j++)
          for (int i = 0; i < C_.rows(); i++)
            adjC(i, j) = C_(i, j).vi_->adj_;
        chainAB(A_, B_, value_of(A_), value_of(B_),
                adjC);
      }

      Eigen::Matrix<TA, RA, CA>  A_;
      Eigen::Matrix<TB, RB, CB>  B_;
      Eigen::Matrix<var, CB, CB> C_;

      Eigen::Map<Eigen::Matrix<TA, RA, CA> > map_A_;
      Eigen::Map<Eigen::Matrix<TB, RB, CB> > map_B_;
      Eigen::Map<Eigen::Matrix<var, CB, CB> > map_C_;

      bool sym_;
    };


    template <typename TA, int RA, int CA, typename TB, int RB, int CB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value
                       || boost::is_same<TB, var>::value,
                       Eigen::Matrix<var, CB, CB> >::type
    quad_form(const Eigen::Matrix<TA, RA, CA>& A,
              const Eigen::Matrix<TB, RB, CB>& B) {
      check_square("quad_form", "A", A);
      check_multiplicable("quad_form", "A", A, "B", B);
      quad_form_vari<TA, RA, CA, TB, RB, CB> *baseVari
        = new quad_form_vari<TA, RA, CA, TB, RB, CB>(A, B);
      return baseVari->C_;
    }

    template <typename TA, int RA, int CA, typename TB, int RB>
    inline typename
    boost::enable_if_c<boost::is_same<TA, var>::value
                       || boost::is_same<TB, var>::value,
                       var>::type
    quad_form(const Eigen::Matrix<TA, RA, CA>& A,
              const Eigen::Matrix<TB, RB, 1>& B) {
      check_square("quad_form", "A", A);
      check_multiplicable("quad_form", "A", A, "B", B);
      quad_form_vari<TA, RA, CA, TB, RB, 1> *baseVari
        = new quad_form_vari<TA, RA, CA, TB, RB, 1>(A, B);
      return baseVari->C_(0, 0);
    }

  }
}
#endif
