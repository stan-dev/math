#ifndef STAN_MATH_REV_MAT_FUN_TRACE_QUAD_FORM_VARI_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_QUAD_FORM_VARI_HPP

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/trace_quad_form.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/rev/mat/fun/adjoints_of.hpp>
#include <stan/math/rev/mat/fun/increment_adjoint.hpp>
#include <stan/math/rev/mat/fun/memalloc_matrix_map.hpp>


namespace stan {
  namespace math {

    template <typename TA, int RA, int CA, typename TB, int RB, int CB>
    class trace_quad_form_vari : public vari {
    private:
      typedef Eigen::Matrix<double, RA, CA> ad_t;
      typedef Eigen::Matrix<double, RB, CB> bd_t;

      typedef Eigen::Matrix<var, RA, CA> av_t;
      typedef Eigen::Matrix<var, RB, CB> bv_t;

      typedef Eigen::Matrix<TA, RA, CA> at_t;
      typedef Eigen::Matrix<TB, RB, CB> bt_t;

    protected:
      static inline void chainA(const Eigen::Map<ad_t>& A,
                                const bd_t& Bd, double adjC) {
      }

      static inline void chainA(const Eigen::Map<av_t>& A, const bd_t& Bd,
                                double adjC) {
        increment_adjoint(A, adjC * Bd * Bd.transpose());
      }

      static inline void chainB(const Eigen::Map<bd_t>& B,
                                const ad_t& Ad, const bd_t& Bd, double adjC) {
      }

      static inline void chainB(const Eigen::Map<bv_t>& B,
                                const ad_t& Ad, const bd_t& Bd, double adjC) {
        increment_adjoint(B, adjC * (Ad + Ad.transpose()) * Bd);
      }

      inline void chainAB(const Eigen::Map<at_t>& A, const Eigen::Map<bt_t>& B,
                          const ad_t& Ad, const bd_t& Bd, double adjC) {
        chainA(A, Bd, adjC);
        chainB(B, Ad, Bd, adjC);
      }

    public:
      trace_quad_form_vari(const Eigen::Matrix<TA, RA, CA>& A,
                           const Eigen::Matrix<TB, RB, CB>& B)
        : vari(trace_quad_form(value_of(A), value_of(B))),
          map_A_(memalloc_matrix_map(A)),
          map_B_(memalloc_matrix_map(B)) { }

      virtual void chain() {
        chainAB(map_A_, map_B_,
                value_of(at_t(map_A_)), value_of(bt_t(map_B_)), adj_);
      }

      Eigen::Map<at_t> map_A_;
      Eigen::Map<bt_t> map_B_;
    };

  }
}
#endif
