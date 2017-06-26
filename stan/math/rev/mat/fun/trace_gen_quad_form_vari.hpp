#ifndef STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_VARI_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_VARI_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/mat/fun/adjoints_of.hpp>
#include <stan/math/rev/mat/fun/increment_adjoint.hpp>
#include <stan/math/rev/mat/fun/memalloc_matrix_map.hpp>

namespace stan {
  namespace math {

    template <typename TD, int RD, int CD, typename TA, int RA, int CA,
              typename TB, int RB, int CB>
    class trace_gen_quad_form_vari : public vari {
    private:
      typedef Eigen::Matrix<double, RA, CA> ad_t;
      typedef Eigen::Matrix<double, RB, CB> bd_t;
      typedef Eigen::Matrix<double, RD, CD> dd_t;

      typedef Eigen::Matrix<var, RA, CA> av_t;
      typedef Eigen::Matrix<var, RB, CB> bv_t;
      typedef Eigen::Matrix<var, RD, CD> dv_t;

      typedef Eigen::Matrix<TA, RA, CA> at_t;
      typedef Eigen::Matrix<TB, RB, CB> bt_t;
      typedef Eigen::Matrix<TD, RB, CB> dt_t;

    protected:
      inline void compute_adjoints_A(const Eigen::Map<ad_t>& Avar,
                                     double adj, const bd_t& B,
                                     const Eigen::Matrix<double, RA, CB>& BD) {
      }

      inline void compute_adjoints_A(const Eigen::Map<av_t>& Avar,
                                     double adj, const bd_t& B,
                                     const Eigen::Matrix<double, RA, CB>& BD) {
        increment_adjoint(Avar, adj * (B * BD.transpose()));
      }

      inline void compute_adjoints_B(const Eigen::Map<bd_t>& Bvar,
                                     double adj, const ad_t& A,
                                     const Eigen::Matrix<double, RA, CB>& BD,
                                     const Eigen::Matrix<double, CA, CB>& AtB,
                                     const ad_t& D) {
      }

      inline void compute_adjoints_B(const Eigen::Map<bv_t>& Bvar,
                                     double adj, const ad_t& A,
                                     const Eigen::Matrix<double, RA, CB>& BD,
                                     const Eigen::Matrix<double, CA, CB>& AtB,
                                     const ad_t& D) {
        increment_adjoint(Bvar, adj * (A * BD + AtB * D.transpose()));
      }

      inline void compute_adjoints_D(const Eigen::Map<dd_t>& Dvar,
                                     double adj, const bd_t& B,
                                     const Eigen::Matrix<double, CA, CB>& AtB,
                                     const ad_t& D) {
      }

      inline void compute_adjoints_D(const Eigen::Map<dv_t>& Dvar,
                                     double adj, const bd_t& B,
                                     const Eigen::Matrix<double, CA, CB>& AtB,
                                     const ad_t& D) {
        increment_adjoint(Dvar, adj * (B.transpose() * AtB));
      }

      inline void compute_adjoints(double adj, const dd_t& D, const ad_t& A,
                                   const bd_t& B) {
        Eigen::Matrix<double, RA, CB> BD;
        if (boost::is_same<TB, var>::value ||  boost::is_same<TA, var>::value)
          BD.noalias() = B * D;

        Eigen::Matrix<double, CA, CB> AtB;
        if (boost::is_same<TB, var>::value ||  boost::is_same<TD, var>::value)
          AtB.noalias() = A.transpose() * B;

        compute_adjoints_B(B_, adj, A, BD, AtB, D);
        compute_adjoints_A(A_, adj, B, BD);
        compute_adjoints_D(D_, adj, B, AtB, D);
      }

    public:
      trace_gen_quad_form_vari(const Eigen::Matrix<TD, RD, CD>& D,
                               const Eigen::Matrix<TA, RA, CA>& A,
                               const Eigen::Matrix<TB, RB, CB>& B)
        : vari(trace_gen_quad_form(value_of(D), value_of(A), value_of(B))),
          D_(memalloc_matrix_map(D)), A_(memalloc_matrix_map(A)),
          B_(memalloc_matrix_map(B)) { }

      virtual void chain() {
        compute_adjoints(adj_, value_of(dt_t(D_)), value_of(at_t(A_)),
                         value_of(bt_t(B_)));
      }

      Eigen::Map<dt_t> D_;
      Eigen::Map<at_t> A_;
      Eigen::Map<bt_t> B_;
    };

  }
}
#endif
