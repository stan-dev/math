#ifndef STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/trace_gen_quad_form_vari.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/fun/trace_gen_quad_form.hpp>

namespace stan {
  namespace math {

    template <typename TD, int RD, int CD, typename TA, int RA, int CA,
              typename TB, int RB, int CB>
    inline typename boost::enable_if_c<boost::is_same<TD, var>::value
                                       || boost::is_same<TA, var>::value
                                       || boost::is_same<TB, var>::value,
                                       var>::type
    trace_gen_quad_form(const Eigen::Matrix<TD, RD, CD>& D,
                        const Eigen::Matrix<TA, RA, CA>& A,
                        const Eigen::Matrix<TB, RB, CB>& B) {
      check_square("trace_gen_quad_form", "A", A);
      check_square("trace_gen_quad_form", "D", D);
      check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
      check_multiplicable("trace_gen_quad_form", "B", B, "D", D);

      return var(new trace_gen_quad_form_vari<TD, RD, CD, TA, RA,
                                              CA, TB, RB, CB>(D, A, B));
    }

}
  }
#endif
