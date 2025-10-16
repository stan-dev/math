#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/mock_ode_functor.hpp>
#include <test/unit/math/prim/functor/mock_throwing_ode_functor.hpp>
#include <vector>
#include <string>

struct StanAgradRevOde : public ::testing::Test {
  void SetUp() { stan::math::recover_memory(); }
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
};

// ******************** DV ****************************
TEST_F(StanAgradRevOde, coupled_ode_system_dv) {
  using stan::math::coupled_ode_system;

  // Run nested autodiff in this scope
  stan::math::nested_rev_autodiff nested;

  harm_osc_ode_fun_eigen harm_osc;

  std::vector<stan::math::var> theta;
  std::vector<double> z0;
  Eigen::VectorXd y0(2);
  double t0;
  const size_t z_size = 4;
  std::vector<double> dz_dt(z_size, 0);

  double gamma(0.15);
  t0 = 0;

  theta.push_back(gamma);
  y0 << 1.0, 0.5;

  z0.push_back(1.0);
  z0.push_back(0.5);
  z0.push_back(1.0);
  z0.push_back(2.0);

  std::size_t stack_size = stan::math::nested_size();

  coupled_ode_system<decltype(harm_osc), double, std::vector<stan::math::var>,
                     std::vector<double>, std::vector<int>>
      system(harm_osc, y0, &msgs, theta, x, x_int);

  EXPECT_EQ(stack_size, stan::math::nested_size())
      << "expecting no new things on the stack";

  system(z0, dz_dt, t0);

  EXPECT_FLOAT_EQ(0.5, dz_dt[0]);
  EXPECT_FLOAT_EQ(-1.075, dz_dt[1]);
  EXPECT_FLOAT_EQ(2, dz_dt[2]);
  EXPECT_FLOAT_EQ(-1.8, dz_dt[3]);
}

TEST_F(StanAgradRevOde, initial_state_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;

  Eigen::VectorXd y0_d = Eigen::VectorXd::Zero(N);
  std::vector<var> theta_v(M, 0.0);

  for (size_t n = 0; n < N; n++)
    y0_d[n] = n + 1;
  for (size_t m = 0; m < M; m++)
    theta_v[m] = 10 * (m + 1);

  coupled_ode_system<decltype(base_ode), double, std::vector<stan::math::var>,
                     std::vector<double>, std::vector<int>>
      coupled_system_dv(base_ode, y0_d, &msgs, theta_v, x, x_int);

  std::vector<double> state = coupled_system_dv.initial_state();
  for (size_t n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(y0_d[n], state[n])
        << "we don't need derivatives of y0; "
        << "initial state gets the initial values";
  for (size_t n = N; n < state.size(); n++)
    EXPECT_FLOAT_EQ(0.0, state[n]);
}

TEST_F(StanAgradRevOde, size_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * M;

  Eigen::VectorXd y0_d = Eigen::VectorXd::Zero(N);
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<decltype(base_ode), double, std::vector<stan::math::var>,
                     std::vector<double>, std::vector<int>>
      coupled_system_dv(base_ode, y0_d, &msgs, theta_v, x, x_int);

  EXPECT_EQ(z_size, coupled_system_dv.size());
}

TEST_F(StanAgradRevOde, memory_recovery_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * M;

  Eigen::VectorXd y0_d = Eigen::VectorXd::Zero(N);
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<decltype(base_ode), double, std::vector<stan::math::var>,
                     std::vector<double>, std::vector<int>>
      coupled_system_dv(base_ode, y0_d, &msgs, theta_v, x, x_int);

  std::vector<double> z(z_size, 0);
  std::vector<double> dz_dt(z_size, 0);
  double t = 10;

  EXPECT_TRUE(stan::math::empty_nested());
  EXPECT_NO_THROW(coupled_system_dv(z, dz_dt, t));
  EXPECT_TRUE(stan::math::empty_nested());
}

TEST_F(StanAgradRevOde, memory_recovery_exception_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  std::string message = "ode throws";

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * M;

  for (size_t n = 0; n < N + 1; n++) {
    std::stringstream scoped_message;
    scoped_message << "iteration " << n;
    SCOPED_TRACE(scoped_message.str());
    mock_throwing_ode_functor<std::logic_error> throwing_ode(message, 1);

    Eigen::VectorXd y0_d = Eigen::VectorXd::Zero(N);
    std::vector<var> theta_v(M, 0.0);

    coupled_ode_system<decltype(throwing_ode), double,
                       std::vector<stan::math::var>, std::vector<double>,
                       std::vector<int>>
        coupled_system_dv(throwing_ode, y0_d, &msgs, theta_v, x, x_int);

    std::vector<double> z(z_size, 0);
    std::vector<double> dz_dt(z_size, 0);
    double t = 10;

    EXPECT_TRUE(stan::math::empty_nested());
    EXPECT_THROW_MSG(coupled_system_dv(z, dz_dt, t), std::logic_error, message);
    EXPECT_TRUE(stan::math::empty_nested());
  }
}

// ******************** VD ****************************

TEST_F(StanAgradRevOde, coupled_ode_system_vd) {
  using stan::math::coupled_ode_system;

  // Run nested autodiff in this scope
  stan::math::nested_rev_autodiff nested;

  harm_osc_ode_fun_eigen harm_osc;

  std::vector<double> theta;
  std::vector<double> z0;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_var(2);
  std::vector<double> y0_adj;
  double t0;
  const size_t N = 2;
  const size_t z_size = N + N * N;
  std::vector<double> dz_dt(z_size, 0);

  double gamma(0.15);
  t0 = 0;

  theta.push_back(gamma);

  z0.push_back(1.0);
  z0.push_back(0.5);
  z0.push_back(1.0);
  z0.push_back(0.0);
  z0.push_back(0.0);
  z0.push_back(1.0);

  y0_var << 1.0, 0.5;

  std::size_t stack_size = stan::math::nested_size();

  coupled_ode_system<decltype(harm_osc), stan::math::var, std::vector<double>,
                     std::vector<double>, std::vector<int>>
      system(harm_osc, y0_var, &msgs, theta, x, x_int);

  EXPECT_EQ(stack_size, stan::math::nested_size())
      << "expecting no new things on the stack";

  system(z0, dz_dt, t0);

  EXPECT_FLOAT_EQ(0.5, dz_dt[0]);
  EXPECT_FLOAT_EQ(-1.0 - 0.15 * 0.5, dz_dt[1]);
  EXPECT_FLOAT_EQ(0.0 * 1.0 + 1.0 * 0.0, dz_dt[2]);
  EXPECT_FLOAT_EQ(-1.0 * 1.0 - 0.15 * 0.0, dz_dt[3]);
  EXPECT_FLOAT_EQ(0.0 * 0.0 + 1.0 * 1.0, dz_dt[4]);
  EXPECT_FLOAT_EQ(-1.0 * 0.0 - 0.15 * 1.0, dz_dt[5]);
}

TEST_F(StanAgradRevOde, initial_state_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
      = Eigen::VectorXd::Zero(N).template cast<var>();
  std::vector<double> theta_d(M, 0.0);

  for (size_t n = 0; n < N; n++)
    y0_v[n] = n + 1;
  for (size_t m = 0; m < M; m++)
    theta_d[m] = 10 * (m + 1);

  coupled_ode_system<decltype(base_ode), stan::math::var, std::vector<double>,
                     std::vector<double>, std::vector<int>>
      coupled_system_vd(base_ode, y0_v, &msgs, theta_d, x, x_int);

  std::vector<double> state;

  state = coupled_system_vd.initial_state();
  for (size_t n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(n + 1, state[n]);
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(i == j ? 1.0 : 0.0, state[N + i + j * N]);
}

TEST_F(StanAgradRevOde, size_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
      = Eigen::VectorXd::Zero(N).template cast<var>();
  std::vector<double> theta_d(M, 0.0);

  coupled_ode_system<decltype(base_ode), stan::math::var, std::vector<double>,
                     std::vector<double>, std::vector<int>>
      coupled_system_vd(base_ode, y0_v, &msgs, theta_d, x, x_int);

  EXPECT_EQ(z_size, coupled_system_vd.size());
}

TEST_F(StanAgradRevOde, memory_recovery_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
      = Eigen::VectorXd::Zero(N).template cast<var>();
  std::vector<double> theta_d(M, 0.0);

  coupled_ode_system<decltype(base_ode), stan::math::var, std::vector<double>,
                     std::vector<double>, std::vector<int>>
      coupled_system_vd(base_ode, y0_v, &msgs, theta_d, x, x_int);

  std::vector<double> z(z_size, 0);
  std::vector<double> dz_dt(z_size, 0);
  double t = 10;

  EXPECT_TRUE(stan::math::empty_nested());
  EXPECT_NO_THROW(coupled_system_vd(z, dz_dt, t));
  EXPECT_TRUE(stan::math::empty_nested());
}

TEST_F(StanAgradRevOde, memory_recovery_exception_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  std::string message = "ode throws";

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N;

  for (size_t n = 0; n < N + 1; n++) {
    std::stringstream scoped_message;
    scoped_message << "iteration " << n;
    SCOPED_TRACE(scoped_message.str());
    mock_throwing_ode_functor<std::logic_error> throwing_ode(message, 1);

    Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
        = Eigen::VectorXd::Zero(N).template cast<var>();
    std::vector<double> theta_d(M, 0.0);

    coupled_ode_system<decltype(throwing_ode), stan::math::var,
                       std::vector<double>, std::vector<double>,
                       std::vector<int>>
        coupled_system_vd(throwing_ode, y0_v, &msgs, theta_d, x, x_int);

    std::vector<double> z(z_size, 0);
    std::vector<double> dz_dt(z_size, 0);
    double t = 10;

    EXPECT_TRUE(stan::math::empty_nested());
    EXPECT_THROW_MSG(coupled_system_vd(z, dz_dt, t), std::logic_error, message);
    EXPECT_TRUE(stan::math::empty_nested());
  }
}

// ******************** VV ****************************

TEST_F(StanAgradRevOde, coupled_ode_system_vv) {
  using stan::math::coupled_ode_system;

  // Run nested autodiff in this scope
  stan::math::nested_rev_autodiff nested;

  const size_t N = 2;
  const size_t M = 1;
  const size_t z_size = N + N * N + N * M;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_var(2);
  y0_var << 1.0, 0.5;

  std::vector<stan::math::var> theta_var(1);
  theta_var[0] = 0.15;

  harm_osc_ode_fun_eigen harm_osc;

  std::size_t stack_size = stan::math::nested_size();

  coupled_ode_system<decltype(harm_osc), stan::math::var,
                     std::vector<stan::math::var>, std::vector<double>,
                     std::vector<int>>
      system(harm_osc, y0_var, &msgs, theta_var, x, x_int);

  EXPECT_EQ(stack_size, stan::math::nested_size())
      << "expecting no new things on the stack";

  std::vector<double> z0(z_size, 0);
  z0[0] = 1.0;
  z0[1] = 0.5;
  z0[2] = 1.0;
  z0[5] = 1.0;

  double t0;
  t0 = 0;

  std::vector<double> dz_dt(z_size);
  system(z0, dz_dt, t0);

  Eigen::VectorXd y0_double(2);
  y0_double << 1.0, 0.5;

  std::vector<double> theta_double(1);
  theta_double[0] = 0.15;

  Eigen::VectorXd dy_dt_base
      = harm_osc(0.0, y0_double, &msgs, theta_double, x, x_int);

  EXPECT_FLOAT_EQ(dy_dt_base[0], dz_dt[0]);
  EXPECT_FLOAT_EQ(dy_dt_base[1], dz_dt[1]);
  EXPECT_FLOAT_EQ(0, dz_dt[2]);
  EXPECT_FLOAT_EQ(-1, dz_dt[3]);
  EXPECT_FLOAT_EQ(1, dz_dt[4]);
  EXPECT_FLOAT_EQ(-0.15, dz_dt[5]);
  EXPECT_FLOAT_EQ(0, dz_dt[6]);
  EXPECT_FLOAT_EQ(-0.5, dz_dt[7]);
}

TEST_F(StanAgradRevOde, initial_state_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
      = Eigen::VectorXd::Zero(N).template cast<var>();
  std::vector<var> theta_v(M, 0.0);

  for (size_t n = 0; n < N; n++)
    y0_v[n] = n + 1;
  for (size_t m = 0; m < M; m++)
    theta_v[m] = 10 * (m + 1);

  coupled_ode_system<decltype(base_ode), stan::math::var,
                     std::vector<stan::math::var>, std::vector<double>,
                     std::vector<int>>
      coupled_system_vv(base_ode, y0_v, &msgs, theta_v, x, x_int);

  std::vector<double> state = coupled_system_vv.initial_state();
  for (size_t n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(n + 1, state[n]);
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(i == j ? 1.0 : 0.0, state[N + i + j * N]);
  for (size_t n = N + N * N; n < N + N * N + N * M; n++)
    EXPECT_FLOAT_EQ(0.0, state[n]);
}

TEST_F(StanAgradRevOde, size_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N + N * M;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
      = Eigen::VectorXd::Zero(N).template cast<var>();
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<decltype(base_ode), stan::math::var,
                     std::vector<stan::math::var>, std::vector<double>,
                     std::vector<int>>
      coupled_system_vv(base_ode, y0_v, &msgs, theta_v, x, x_int);

  EXPECT_EQ(z_size, coupled_system_vv.size());
}

TEST_F(StanAgradRevOde, memory_recovery_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N + N * M;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
      = Eigen::VectorXd::Zero(N).template cast<var>();
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<decltype(base_ode), stan::math::var,
                     std::vector<stan::math::var>, std::vector<double>,
                     std::vector<int>>
      coupled_system_vv(base_ode, y0_v, &msgs, theta_v, x, x_int);

  std::vector<double> z(z_size, 0);
  std::vector<double> dz_dt(z_size, 0);
  double t = 10;

  EXPECT_TRUE(stan::math::empty_nested());
  EXPECT_NO_THROW(coupled_system_vv(z, dz_dt, t));
  EXPECT_TRUE(stan::math::empty_nested());
}

TEST_F(StanAgradRevOde, memory_recovery_exception_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  std::string message = "ode throws";

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N + N * M;

  for (size_t n = 0; n < N + 1; n++) {
    std::stringstream scoped_message;
    scoped_message << "iteration " << n;
    SCOPED_TRACE(scoped_message.str());
    mock_throwing_ode_functor<std::logic_error> throwing_ode(message, 1);

    Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y0_v
        = Eigen::VectorXd::Zero(N).template cast<var>();
    std::vector<var> theta_v(M, 0.0);

    coupled_ode_system<decltype(throwing_ode), stan::math::var,
                       std::vector<stan::math::var>, std::vector<double>,
                       std::vector<int>>
        coupled_system_vv(throwing_ode, y0_v, &msgs, theta_v, x, x_int);

    std::vector<double> z(z_size, 0);
    std::vector<double> dz_dt(z_size, 0);
    double t = 10;

    EXPECT_TRUE(stan::math::empty_nested());
    EXPECT_THROW_MSG(coupled_system_vv(z, dz_dt, t), std::logic_error, message);
    EXPECT_TRUE(stan::math::empty_nested());
  }
}

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

TEST_F(StanAgradRevOde, coupled_ode_system_var) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  double a = 1.3;
  var av = a;

  ayt func;

  coupled_ode_system<ayt, stan::math::var, decltype(av)> system(func, y0v,
                                                                &msgs, av);

  std::vector<double> z = {y0(0), y0(1), 3.2, 3.3, 4.4, 4.5, 2.1, 2.2};

  double t0 = 0.75;

  std::vector<double> dz_dt;
  system(z, dz_dt, t0);

  EXPECT_FLOAT_EQ(dz_dt[0], -0.0975);
  EXPECT_FLOAT_EQ(dz_dt[1], -0.3315);
  EXPECT_FLOAT_EQ(dz_dt[2], -3.12000);
  EXPECT_FLOAT_EQ(dz_dt[3], -5.46975);
  EXPECT_FLOAT_EQ(dz_dt[4], -4.29000);
  EXPECT_FLOAT_EQ(dz_dt[5], -7.45875);
  EXPECT_FLOAT_EQ(dz_dt[6], -2.1225);
  EXPECT_FLOAT_EQ(dz_dt[7], -3.9015);
}

TEST_F(StanAgradRevOde, coupled_ode_system_std_vector) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  std::vector<var> av = {1.3};

  ayt func;

  coupled_ode_system<ayt, stan::math::var, decltype(av)> system(func, y0v,
                                                                &msgs, av);

  std::vector<double> z = {y0(0), y0(1), 3.2, 3.3, 4.4, 4.5, 2.1, 2.2};

  double t0 = 0.75;

  std::vector<double> dz_dt;
  system(z, dz_dt, t0);

  EXPECT_FLOAT_EQ(dz_dt[0], -0.0975);
  EXPECT_FLOAT_EQ(dz_dt[1], -0.3315);
  EXPECT_FLOAT_EQ(dz_dt[2], -3.12000);
  EXPECT_FLOAT_EQ(dz_dt[3], -5.46975);
  EXPECT_FLOAT_EQ(dz_dt[4], -4.29000);
  EXPECT_FLOAT_EQ(dz_dt[5], -7.45875);
  EXPECT_FLOAT_EQ(dz_dt[6], -2.1225);
  EXPECT_FLOAT_EQ(dz_dt[7], -3.9015);
}

TEST_F(StanAgradRevOde, coupled_ode_system_vector) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  Eigen::Matrix<var, Eigen::Dynamic, 1> av(1);
  av << 1.3;

  ayt func;

  coupled_ode_system<ayt, stan::math::var, decltype(av)> system(func, y0v,
                                                                &msgs, av);

  std::vector<double> z = {y0(0), y0(1), 3.2, 3.3, 4.4, 4.5, 2.1, 2.2};

  double t0 = 0.75;

  std::vector<double> dz_dt;
  system(z, dz_dt, t0);

  EXPECT_FLOAT_EQ(dz_dt[0], -0.0975);
  EXPECT_FLOAT_EQ(dz_dt[1], -0.3315);
  EXPECT_FLOAT_EQ(dz_dt[2], -3.12000);
  EXPECT_FLOAT_EQ(dz_dt[3], -5.46975);
  EXPECT_FLOAT_EQ(dz_dt[4], -4.29000);
  EXPECT_FLOAT_EQ(dz_dt[5], -7.45875);
  EXPECT_FLOAT_EQ(dz_dt[6], -2.1225);
  EXPECT_FLOAT_EQ(dz_dt[7], -3.9015);
}

TEST_F(StanAgradRevOde, coupled_ode_system_row_vector) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  Eigen::Matrix<var, 1, Eigen::Dynamic> av(1);
  av << 1.3;

  ayt func;

  coupled_ode_system<ayt, stan::math::var, decltype(av)> system(func, y0v,
                                                                &msgs, av);

  std::vector<double> z = {y0(0), y0(1), 3.2, 3.3, 4.4, 4.5, 2.1, 2.2};

  double t0 = 0.75;

  std::vector<double> dz_dt;
  system(z, dz_dt, t0);

  EXPECT_FLOAT_EQ(dz_dt[0], -0.0975);
  EXPECT_FLOAT_EQ(dz_dt[1], -0.3315);
  EXPECT_FLOAT_EQ(dz_dt[2], -3.12000);
  EXPECT_FLOAT_EQ(dz_dt[3], -5.46975);
  EXPECT_FLOAT_EQ(dz_dt[4], -4.29000);
  EXPECT_FLOAT_EQ(dz_dt[5], -7.45875);
  EXPECT_FLOAT_EQ(dz_dt[6], -2.1225);
  EXPECT_FLOAT_EQ(dz_dt[7], -3.9015);
}

TEST_F(StanAgradRevOde, coupled_ode_system_matrix) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> av(1, 1);
  av << 1.3;

  ayt func;

  coupled_ode_system<ayt, stan::math::var, decltype(av)> system(func, y0v,
                                                                &msgs, av);

  std::vector<double> z = {y0(0), y0(1), 3.2, 3.3, 4.4, 4.5, 2.1, 2.2};

  double t0 = 0.75;

  std::vector<double> dz_dt;
  system(z, dz_dt, t0);

  EXPECT_FLOAT_EQ(dz_dt[0], -0.0975);
  EXPECT_FLOAT_EQ(dz_dt[1], -0.3315);
  EXPECT_FLOAT_EQ(dz_dt[2], -3.12000);
  EXPECT_FLOAT_EQ(dz_dt[3], -5.46975);
  EXPECT_FLOAT_EQ(dz_dt[4], -4.29000);
  EXPECT_FLOAT_EQ(dz_dt[5], -7.45875);
  EXPECT_FLOAT_EQ(dz_dt[6], -2.1225);
  EXPECT_FLOAT_EQ(dz_dt[7], -3.9015);
}

TEST_F(StanAgradRevOde, coupled_ode_system_extra_args) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  Eigen::VectorXd y0(2);
  y0 << 0.1, 0.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0v = y0.template cast<var>();

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> av(1, 1);
  av << -0.2;
  int e1 = 1;
  double e2 = 0.1;
  std::vector<double> e3 = {0.1};
  Eigen::VectorXd e4(1);
  e4 << 0.1;
  Eigen::VectorXd e5(1);
  e5 << 0.1;
  Eigen::MatrixXd e6(1, 1);
  e6 << 0.1;

  ayt func;

  coupled_ode_system<ayt, stan::math::var, decltype(av), decltype(e1),
                     decltype(e2), decltype(e3), decltype(e4), decltype(e5),
                     decltype(e6)>
      system(func, y0v, &msgs, av, e1, e2, e3, e4, e5, e6);

  std::vector<double> z = {y0(0), y0(1), 3.2, 3.3, 4.4, 4.5, 2.1, 2.2};

  double t0 = 0.75;

  std::vector<double> dz_dt;
  system(z, dz_dt, t0);

  EXPECT_FLOAT_EQ(dz_dt[0], -0.0975);
  EXPECT_FLOAT_EQ(dz_dt[1], -0.3315);
  EXPECT_FLOAT_EQ(dz_dt[2], -3.12000);
  EXPECT_FLOAT_EQ(dz_dt[3], -5.46975);
  EXPECT_FLOAT_EQ(dz_dt[4], -4.29000);
  EXPECT_FLOAT_EQ(dz_dt[5], -7.45875);
  EXPECT_FLOAT_EQ(dz_dt[6], -2.1225);
  EXPECT_FLOAT_EQ(dz_dt[7], -3.9015);
}
