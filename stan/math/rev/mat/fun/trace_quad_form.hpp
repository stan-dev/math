#ifndef STAN_MATH_REV_MAT_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_QUAD_FORM_HPP

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/rev/mat/fun/trace_quad_form_vari.hpp>
#include <stan/math/prim/mat/fun/trace_quad_form.hpp>

namespace stan {
  namespace math {

    template <typename TA, int RA, int CA, typename TB, int RB, int CB>
    inline typename boost::enable_if_c<boost::is_same<TA, var>::value
                                       || boost::is_same<TB, var>::value,
                                       var>::type
    trace_quad_form(const Eigen::Matrix<TA, RA, CA>& A,
                    const Eigen::Matrix<TB, RB, CB>& B) {
      check_square("trace_quad_form", "A", A);
      check_multiplicable("trace_quad_form", "A", A, "B", B);
      return var(new trace_quad_form_vari<TA, RA, CA, TB, RB, CB>(A, B));
    }

  }
}
#endif
