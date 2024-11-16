#ifndef TEST_UNIT_MATH_REV_FUN_ONE_CPT_ODE_HPP
#define TEST_UNIT_MATH_REV_FUN_ONE_CPT_ODE_HPP

#include <stan/math/rev.hpp>
#include <vector>

namespace stan {

namespace math {

  struct one_cpt_ode {
    template <typename T, typename T1, typename T2, typename Tr>
    inline Eigen::Matrix<typename stan::return_type_t<T, T1, T2, Tr>, -1, 1>
    operator()(const Eigen::Matrix<T, -1, 1>& y,
	       const T1& ka, const T2& k10,
	       const Eigen::Matrix<Tr, -1, 1>&) const {
      typedef typename stan::return_type<T, T1, T2, Tr>::type scalar;
      Eigen::Matrix<scalar, -1, 1> yd(2);
      yd << -ka * y(0), ka * y(0) - k10 * y(1);

      return y;
    }
  };


}  // namespace math
}  // namespace stan
#endif
