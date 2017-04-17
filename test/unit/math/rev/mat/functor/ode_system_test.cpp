#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>

struct StanMathRevOdeSystem : public ::testing::Test {
  void SetUp() {
    t0 = 0;

    y0.push_back(1);
    y0.push_back(2);

    theta.push_back(0.5);

    N = 2;
    M = 1;

    Jy_ref.resize(N,N);
    Jy_ref(0,0) =  0;
    Jy_ref(0,1) =  1;
    Jy_ref(1,0) = -1;
    Jy_ref(1,1) = -theta[0];

    Jtheta_ref.resize(N,M);
    Jtheta_ref(0,0) = 0;
    Jtheta_ref(1,0) = -y0[1];
  }
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
  double t0;
  std::vector<double> y0;
  std::vector<double> theta;
  harm_osc_ode_fun ode_rhs;
  size_t N;
  size_t M;
  Eigen::MatrixXd Jy_ref;
  Eigen::MatrixXd Jtheta_ref;
};


// ************ ODE RHS functor ************************

TEST_F(StanMathRevOdeSystem, ode_system_rhs) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  std::vector<double> dy_dt;

  ode_system(t0, y0, dy_dt);

  EXPECT_FLOAT_EQ( 2.0, dy_dt[0]);
  EXPECT_FLOAT_EQ(-2.0, dy_dt[1]);
}

// ************ jacobian wrt to states ************************

TEST_F(StanMathRevOdeSystem, ode_system_jac_v_m) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  Eigen::VectorXd dy_dt(N);
  Eigen::MatrixXd Jy(N,N);

  ode_system.jacobian(t0, y0, dy_dt, Jy);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

}

TEST_F(StanMathRevOdeSystem, ode_system_jac_Mv_m) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  std::vector<double> dy_dt_raw(N, 0);
  Eigen::Map<Eigen::VectorXd> dy_dt(&dy_dt_raw[0], N);
  Eigen::MatrixXd Jy(N,N);

  ode_system.jacobian(t0, y0, dy_dt, Jy);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

}

TEST_F(StanMathRevOdeSystem, ode_system_jac_v_Mm) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  Eigen::VectorXd dy_dt(N);
  std::vector<double> Jy_raw(N*N, 0);
  Eigen::Map<Eigen::MatrixXd> Jy(&Jy_raw[0],N,N);

  ode_system.jacobian(t0, y0, dy_dt, Jy);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

}

TEST_F(StanMathRevOdeSystem, ode_system_jac_Mv_Mm) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  std::vector<double> dy_dt_raw(N, 0);
  Eigen::Map<Eigen::VectorXd> dy_dt(&dy_dt_raw[0], N);
  std::vector<double> Jy_raw(N*N, 0);
  Eigen::Map<Eigen::MatrixXd> Jy(&Jy_raw[0],N,N);

  ode_system.jacobian(t0, y0, dy_dt, Jy);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

}


// ********** jacobian wrt to states + parameters *************

TEST_F(StanMathRevOdeSystem, ode_system_jac_v_m_m) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  Eigen::VectorXd dy_dt(N);
  Eigen::MatrixXd Jy(N,N);
  Eigen::MatrixXd Jtheta(N,M);

  ode_system.jacobian(t0, y0, dy_dt, Jy, Jtheta);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < M; j++)
      EXPECT_FLOAT_EQ(Jtheta_ref(i,j), Jtheta(i,j));
}

TEST_F(StanMathRevOdeSystem, ode_system_jac_Mv_m_m) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  std::vector<double> dy_dt_raw(N, 0);
  Eigen::Map<Eigen::VectorXd> dy_dt(&dy_dt_raw[0], N);
  Eigen::MatrixXd Jy(N,N);
  Eigen::MatrixXd Jtheta(N,M);

  ode_system.jacobian(t0, y0, dy_dt, Jy, Jtheta);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < M; j++)
      EXPECT_FLOAT_EQ(Jtheta_ref(i,j), Jtheta(i,j));
}

TEST_F(StanMathRevOdeSystem, ode_system_jac_v_Mm_Mm) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  Eigen::VectorXd dy_dt(N);
  std::vector<double> Jy_raw(N*N, 0);
  Eigen::Map<Eigen::MatrixXd> Jy(&Jy_raw[0],N,N);
  std::vector<double> Jtheta_raw(N*M, 0);
  Eigen::Map<Eigen::MatrixXd> Jtheta(&Jtheta_raw[0],N,M);

  ode_system.jacobian(t0, y0, dy_dt, Jy, Jtheta);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < M; j++)
      EXPECT_FLOAT_EQ(Jtheta_ref(i,j), Jtheta(i,j));
}


TEST_F(StanMathRevOdeSystem, ode_system_jac_Mv_Mm_Mm) {
  stan::math::ode_system<harm_osc_ode_fun> ode_system(ode_rhs, theta, x, x_int, &msgs);

  std::vector<double> dy_dt_raw(N, 0);
  Eigen::Map<Eigen::VectorXd> dy_dt(&dy_dt_raw[0], N);
  std::vector<double> Jy_raw(N*N, 0);
  Eigen::Map<Eigen::MatrixXd> Jy(&Jy_raw[0],N,N);
  std::vector<double> Jtheta_raw(N*M, 0);
  Eigen::Map<Eigen::MatrixXd> Jtheta(&Jtheta_raw[0],N,M);

  ode_system.jacobian(t0, y0, dy_dt, Jy, Jtheta);

  EXPECT_FLOAT_EQ( 2.0, dy_dt(0));
  EXPECT_FLOAT_EQ(-2.0, dy_dt(1));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(Jy_ref(i,j), Jy(i,j));

  for(size_t i = 0; i < N; i++)
    for(size_t j = 0; j < M; j++)
      EXPECT_FLOAT_EQ(Jtheta_ref(i,j), Jtheta(i,j));
}

