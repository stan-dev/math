#ifndef STAN_MATH_REV_MAT_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_REV_MAT_FUN_QUAD_FORM_SYM_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/rev/mat/fun/quad_form_vari.hpp>

namespace stan {
  namespace math {

    template <typename TA, int RA, int CA, typename TB, int RB, int CB>
    inline typename boost::enable_if_c<boost::is_same<TA, var>::value
                                       || boost::is_same<TB, var>::value,
                                       Eigen::Matrix<var, CB, CB> >::type
    quad_form_sym(const Eigen::Matrix<TA, RA, CA>& A,
                  const Eigen::Matrix<TB, RB, CB>& B) {
      check_square("quad_form_sym", "A", A);
      check_symmetric("quad_form_sym", "A", A);
      check_multiplicable("quad_form_sym", "A", A, "B", B);

      typedef quad_form_vari<TA, RA, CA, TB, RB, CB, true> vari_t;
      vari_t* baseVari = new vari_t(A, B);
      return Eigen::Matrix<var, CB, CB>(baseVari->map_C_);
    }

    template <typename TA, int RA, int CA, typename TB, int RB>
    inline typename boost::enable_if_c<boost::is_same<TA, var>::value
                                       || boost::is_same<TB, var>::value,
                                        var>::type
    quad_form_sym(const Eigen::Matrix<TA, RA, CA>& A,
                  const Eigen::Matrix<TB, RB, 1>& B) {
      check_square("quad_form_sym", "A", A);
      check_symmetric("quad_form_sym", "A", A);
      check_multiplicable("quad_form_sym", "A", A, "B", B);

      typedef quad_form_vari<TA, RA, CA, TB, RB, 1, true> vari_t;
      vari_t* base_vari = new vari_t(A, B);
      return base_vari->map_C_(0, 0);
    }

  }
}
#endif
