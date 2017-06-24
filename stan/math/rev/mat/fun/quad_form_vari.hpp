#ifndef STAN_MATH_REV_MAT_FUN_QUAD_FORM_VARI_HPP
#define STAN_MATH_REV_MAT_FUN_QUAD_FORM_VARI_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>

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

    template <typename TA, int RA, int CA, typename TB, int RB, int CB,
              bool Sym>
    class quad_form_vari : public vari {
    private:
      typedef Eigen::Matrix<double, RA, CA> ad_t;
      typedef Eigen::Matrix<double, RB, CB> bd_t;
      typedef Eigen::Matrix<double, CB, CB> cd_t;

      typedef Eigen::Matrix<var, RA, CA> av_t;
      typedef Eigen::Matrix<var, RB, CB> bv_t;
      typedef Eigen::Matrix<var, CB, CB> cv_t;

    protected:
      inline void chainA(const Eigen::Map<ad_t>& A,
                         const bd_t& Bd, const cd_t& adjC) {
      }

      inline void chainA(const Eigen::Map<av_t>& A,
                         const bd_t& Bd, const cd_t& adjC) {
        Eigen::Matrix<double, RA, CA> adjA(Bd * adjC * Bd.transpose());
        for (int j = 0; j < A.cols(); j++) {
          for (int i = 0; i < A.rows(); i++) {
            A(i, j).vi_->adj_ += adjA(i, j);
          }
        }
      }

      inline void chainB(const Eigen::Map<bd_t>& B,
                         const ad_t& Ad, const bd_t& Bd, const cd_t& adjC) {
      }

      inline void chainB(const Eigen::Map<bv_t>& B,
                         const ad_t& Ad, const bd_t& Bd, const cd_t& adjC) {
        ad_t adjB(Ad * Bd * adjC.transpose() + Ad.transpose() * Bd * adjC);
        for (int j = 0; j < B.cols(); j++)
          for (int i = 0; i < B.rows(); i++)
            B(i, j).vi_->adj_ += adjB(i, j);
      }

      inline void chainAB(const Eigen::Map<Eigen::Matrix<TA, RA, CA> >& A,
                          const Eigen::Map<Eigen::Matrix<TB, RB, CB> >& B,
                          const ad_t& Ad, const bd_t& Bd, const cd_t& adjC) {
        chainA(A, Bd, adjC);
        chainB(B, Ad, Bd, adjC);
      }

    public:
      quad_form_vari(const Eigen::Matrix<TA, RA, CA>& A,
                     const Eigen::Matrix<TB, RB, CB>& B)
        : vari(0.0),
          map_A_(memalloc_matrix_map(A)),
          map_B_(memalloc_matrix_map(B)),
          map_C_(ChainableStack::memalloc_
                   .template alloc_array<var>(B.cols() * B.cols()),
                 B.cols(), B.cols()) {
        bd_t Bd(value_of(B));
        cd_t Cd(Bd.transpose() * value_of(A) * Bd);
        if (Sym) {
          for (int j = 0; j < map_C_.cols(); j++)
            for (int i = 0; i < map_C_.rows(); i++)
              map_C_(i, j) = var(new vari(0.5 * (Cd(i, j) + Cd(j, i)), false));
        } else {
          for (int j = 0; j < map_C_.cols(); j++)
            for (int i = 0; i < map_C_.rows(); i++)
              map_C_(i, j) = var(new vari(Cd(i, j), false));
        }
      }

      virtual void chain() {
        Eigen::Matrix<double, CB, CB> adjC(map_C_.rows(), map_C_.cols());
        for (int j = 0; j < map_C_.cols(); j++)
          for (int i = 0; i < map_C_.rows(); i++)
            adjC(i, j) = map_C_(i, j).vi_->adj_;
        chainAB(map_A_, map_B_,
                value_of(Eigen::Matrix<TA, RA, CA>(map_A_)),
                value_of(Eigen::Matrix<TB, RB, CB>(map_B_)),
                adjC);
      }

      Eigen::Map<Eigen::Matrix<TA, RA, CA> > map_A_;
      Eigen::Map<Eigen::Matrix<TB, RB, CB> > map_B_;
      Eigen::Map<Eigen::Matrix<var, CB, CB> > map_C_;
    };

  }
}
#endif
