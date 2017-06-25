#ifndef STAN_MATH_REV_MAT_FUN_QUAD_FORM_VARI_HPP
#define STAN_MATH_REV_MAT_FUN_QUAD_FORM_VARI_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/adjoints_of.hpp>
#include <stan/math/rev/mat/fun/increment_adjoint.hpp>
#include <stan/math/rev/mat/fun/memalloc_matrix_map.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>

namespace stan {
  namespace math {

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

      typedef Eigen::Matrix<TA, RA, CA> at_t;
      typedef Eigen::Matrix<TB, RB, CB> bt_t;

    protected:
      inline void chainA(const Eigen::Map<ad_t>& A,
                         const bd_t& Bd, const cd_t& adjC) {
      }

      inline void chainA(const Eigen::Map<av_t>& A,
                         const bd_t& Bd, const cd_t& adjC) {
        increment_adjoint(A, Bd * adjC * Bd.transpose());
      }

      inline void chainB(const Eigen::Map<bd_t>& B,
                         const ad_t& Ad, const bd_t& Bd, const cd_t& adjC) {
      }

      inline void chainB(const Eigen::Map<bv_t>& B,
                         const ad_t& Ad, const bd_t& Bd, const cd_t& adjC) {
        increment_adjoint(B, Ad * Bd * adjC.transpose()
                             + Ad.transpose() * Bd * adjC);
      }

      inline void chainAB(const Eigen::Map<at_t>& A, const Eigen::Map<bt_t>& B,
                          const ad_t& Ad, const bd_t& Bd, const cd_t& adjC) {
        chainA(A, Bd, adjC);
        chainB(B, Ad, Bd, adjC);
      }

    public:
      quad_form_vari(const at_t& A, const bt_t& B)
        : vari(0.0),
          map_A_(memalloc_matrix_map(A)),
          map_B_(memalloc_matrix_map(B)),
          map_C_(ChainableStack::memalloc_
                   .template alloc_array<var>(B.cols() * B.cols()),
                 B.cols(), B.cols()) {
        bd_t Bd(value_of(B));
        cd_t Cd(Bd.transpose() * value_of(A) * Bd);
        for (int j = 0; j < map_C_.cols(); ++j)
          for (int i = 0; i < map_C_.rows(); ++i) {
            if (Sym)
              map_C_(i, j) = var(new vari(0.5 * (Cd(i, j) + Cd(j, i)), false));
            else
              map_C_(i, j) = var(new vari(Cd(i, j), false));
          }
      }

      virtual void chain() {
        chainAB(map_A_, map_B_,
                value_of(at_t(map_A_)), value_of(bt_t(map_B_)),
                adjoints_of<CB, CB>(map_C_));
      }

      Eigen::Map<at_t> map_A_;
      Eigen::Map<bt_t> map_B_;
      Eigen::Map<cv_t> map_C_;
    };

  }
}
#endif
