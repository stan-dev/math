#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>
#include <string>

template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
T sum_(T arg) {
  return arg;
}

template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
auto sum_(EigMat&& arg) {
  return stan::math::sum(arg);
}

template <typename Vec, stan::require_std_vector_t<Vec>* = nullptr>
auto sum_(Vec&& arg) {
  stan::scalar_type_t<Vec> sum = 0;
  for (size_t i = 0; i < arg.size(); ++i) {
    sum += sum_(arg[i]);
  }
  return sum;
}

struct ayt {
  template <typename T0, typename T_y, typename... T_Args>
  inline auto operator()(const T0& t, const T_y& y, std::ostream* msgs,
                         const T_Args&... args) const {
    std::vector<typename stan::return_type<T_Args...>::type> vec
        = {sum_(args)...};
    Eigen::Matrix<stan::return_type_t<T0, T_y, T_Args...>, Eigen::Dynamic, 1>
        out(2);
    out(0) = -sum_(vec) * y(0) * t;
    out(1) = -1.7 * sum_(vec) * y(1) * t;
    return out;
  }
};

struct aytm {
  template <typename T0, typename T_y, typename T_Arg, int R, int C>
  inline auto operator()(const T0& t, const T_y& y, std::ostream* msgs,
                         const Eigen::Matrix<T_Arg, R, C>& args) const {
    std::vector<typename stan::return_type<T_Arg>::type> vec = {sum_(args)};
    Eigen::Matrix<stan::return_type_t<T0, T_y, T_Arg>, Eigen::Dynamic, 1> out(
        2);
    out(0) = -sum_(vec) * y(0) * t;
    out(1) = -1.7 * sum_(vec) * y(1) * t;
    return out;
  }
};

TEST(AgradRev, ode_store_sensitivities) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  double a = 1.3;
  var av = a;

  ayt func;

  double t0 = 0.5;
  double t = 0.75;

  var t0v = t0;
  var tv = t;

  std::vector<double> coupled_state = {-0.0975,  -0.3315,  -3.12000, -5.46975,
                                       -4.29000, -7.45875, -2.1225,  -3.9015};

  auto output = stan::math::ode_store_sensitivities(func, coupled_state, y0v,
                                                    t0v, tv, nullptr, av);

  output(0).grad();

  EXPECT_FLOAT_EQ(output(0).val(), -0.0975);
  EXPECT_FLOAT_EQ(y0v(0).adj(), -3.12);
  EXPECT_FLOAT_EQ(y0v(1).adj(), -4.29);
  EXPECT_FLOAT_EQ(av.adj(), -2.1225);
  EXPECT_FLOAT_EQ(t0v.adj(), -1.15089);
  EXPECT_FLOAT_EQ(tv.adj(), 0.0950625);

  stan::math::set_zero_all_adjoints();

  output(1).grad();

  EXPECT_FLOAT_EQ(output(1).val(), -0.3315);
  EXPECT_FLOAT_EQ(y0v(0).adj(), -5.46975);
  EXPECT_FLOAT_EQ(y0v(1).adj(), -7.45875);
  EXPECT_FLOAT_EQ(av.adj(), -3.9015);
  EXPECT_FLOAT_EQ(t0v.adj(), -2.003918);
  EXPECT_FLOAT_EQ(tv.adj(), 0.5494613);

  stan::math::recover_memory();
}

TEST(AgradRev, ode_store_sensitivities_matrix) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  double a = 1.3;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> av(1, 1);
  av(0, 0) = a;

  aytm func;

  double t0 = 0.5;
  double t = 0.75;

  var t0v = t0;
  var tv = t;

  std::vector<double> coupled_state = {-0.0975,  -0.3315,  -3.12000, -5.46975,
                                       -4.29000, -7.45875, -2.1225,  -3.9015};

  auto output = stan::math::ode_store_sensitivities(func, coupled_state, y0v,
                                                    t0v, tv, nullptr, av);

  output(0).grad();

  EXPECT_FLOAT_EQ(output(0).val(), -0.0975);
  EXPECT_FLOAT_EQ(y0v(0).adj(), -3.12);
  EXPECT_FLOAT_EQ(y0v(1).adj(), -4.29);
  EXPECT_FLOAT_EQ(av(0, 0).adj(), -2.1225);
  EXPECT_FLOAT_EQ(t0v.adj(), -1.15089);
  EXPECT_FLOAT_EQ(tv.adj(), 0.0950625);

  stan::math::set_zero_all_adjoints();

  output(1).grad();

  EXPECT_FLOAT_EQ(output(1).val(), -0.3315);
  EXPECT_FLOAT_EQ(y0v(0).adj(), -5.46975);
  EXPECT_FLOAT_EQ(y0v(1).adj(), -7.45875);
  EXPECT_FLOAT_EQ(av(0, 0).adj(), -3.9015);
  EXPECT_FLOAT_EQ(t0v.adj(), -2.003918);
  EXPECT_FLOAT_EQ(tv.adj(), 0.5494613);

  stan::math::recover_memory();
}
