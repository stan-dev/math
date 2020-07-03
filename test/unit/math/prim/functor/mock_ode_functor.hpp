#ifndef TEST_UNIT_MATH_PRIM_FUNCTOR_MOCK_ODE_FUNCTOR
#define TEST_UNIT_MATH_PRIM_FUNCTOR_MOCK_ODE_FUNCTOR
#include <vector>

struct mock_ode_functor {
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2, T3, T4>, Eigen::Dynamic, 1>
  operator()(const T0& t_in, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y_in,
             std::ostream* msgs, const T2& a1, const T3& a2,
             const T4& a3) const {
    return y_in.template cast<stan::return_type_t<T1, T2, T3, T4>>();
  }
};
#endif
