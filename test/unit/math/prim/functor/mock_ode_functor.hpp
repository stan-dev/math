#ifndef TEST_UNIT_MATH_PRIM_FUNCTOR_MOCK_ODE_FUNCTOR
#define TEST_UNIT_MATH_PRIM_FUNCTOR_MOCK_ODE_FUNCTOR
#include <vector>

struct mock_ode_functor {
  template <typename T0, typename T1, typename T2>
  inline std::vector<stan::return_type_t<T1, T2>> operator()(
      const T0& t_in, const std::vector<T1>& y_in, const std::vector<T2>& theta,
      const std::vector<double>& x, const std::vector<int>& x_int,
      std::ostream* msgs) const {
    std::vector<stan::return_type_t<T1, T2>> states;
    states = y_in;
    return states;
  }
};
#endif
