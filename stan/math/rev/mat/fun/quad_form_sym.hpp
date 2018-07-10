#ifndef STAN_MATH_REV_MAT_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_REV_MAT_FUN_QUAD_FORM_SYM_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/quad_form.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/rev/mat/fun/quad_form.hpp>

namespace stan {
namespace math {

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
inline typename boost::enable_if_c<boost::is_same<Ta, var>::value
                                       || boost::is_same<Tb, var>::value,
                                   Eigen::Matrix<var, Cb, Cb> >::type
quad_form_sym(const Eigen::Matrix<Ta, Ra, Ca>& A,
              const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("quad_form", "A", A);
  check_symmetric("quad_form_sym", "A", A);
  check_multiplicable("quad_form_sym", "A", A, "B", B);

  quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>* baseVari
      = new quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>(A, B, true);

  Eigen::Matrix<var, Cb, Cb> C(baseVari->C_rows_, baseVari->C_rows_);
  for (int i = 0; i < baseVari->C_rows_ * baseVari->C_rows_; ++i)
    C.data()[i] = baseVari->C_mem_[i];

  return C;
}

template <typename Ta, int Ra, int Ca, typename Tb, int Rb>
inline typename boost::enable_if_c<
    boost::is_same<Ta, var>::value || boost::is_same<Tb, var>::value, var>::type
quad_form_sym(const Eigen::Matrix<Ta, Ra, Ca>& A,
              const Eigen::Matrix<Tb, Rb, 1>& B) {
  check_square("quad_form", "A", A);
  check_symmetric("quad_form_sym", "A", A);
  check_multiplicable("quad_form_sym", "A", A, "B", B);

  quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>* baseVari
      = new quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>(A, B, true);

  return baseVari->C_mem_[0];
}

}  // namespace math
}  // namespace stan
#endif
