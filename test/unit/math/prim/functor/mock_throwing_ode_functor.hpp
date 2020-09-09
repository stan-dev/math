#ifndef TEST_UNIT_MATH_PRIM_FUNCTOR_MOCK_THROWING_ODE_FUNCTOR
#define TEST_UNIT_MATH_PRIM_FUNCTOR_MOCK_THROWING_ODE_FUNCTOR
#include <vector>
#include <string>

int mock_throwing_ode_functor_count = 0;

template <typename E>
struct mock_throwing_ode_functor {
  const std::string msg_;
  const int N_;  // throw on the N_th call

  explicit mock_throwing_ode_functor(std::string msg) : msg_(msg), N_(1) {
    mock_throwing_ode_functor_count = 0;
  }

  mock_throwing_ode_functor(std::string msg, int N) : msg_(msg), N_(N) {
    mock_throwing_ode_functor_count = 0;
  }

  template <typename T0, typename T1, typename T2>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T0& t_in, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y_in,
             std::ostream* msgs, const std::vector<T2>& theta,
             const std::vector<double>& x,
             const std::vector<int>& x_int) const {
    mock_throwing_ode_functor_count++;
    if (N_ == mock_throwing_ode_functor_count)
      throw E(msg_);
    return y_in;
  }
};
#endif
