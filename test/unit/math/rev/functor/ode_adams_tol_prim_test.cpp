#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/ode_test_functors.hpp>
#include <iostream>
#include <vector>

TEST(ode_adams_tol_prim, y0_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd y0inf(1);
  Eigen::VectorXd y0NaN(1);
  Eigen::VectorXd y0_empty;
  y0NaN << stan::math::NOT_A_NUMBER;
  y0inf << stan::math::INFTY;
  int t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0inf, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0NaN, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0_empty, t0,
                                         ts, 1e-10, 1e-10, 1e6, nullptr, a),
               std::invalid_argument);
}

TEST(ode_adams_tol_prim, t0_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  double t0inf = stan::math::INFTY;
  double t0NaN = stan::math::NOT_A_NUMBER;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0inf, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0NaN, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a),
               std::domain_error);
}

TEST(ode_adams_tol_prim, ts_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};
  std::vector<double> ts_repeat = {0.45, 0.45};
  std::vector<double> ts_lots = {0.45, 0.45, 1.1, 1.1, 2.0};
  std::vector<double> ts_empty = {};
  std::vector<double> ts_early = {-0.45, 0.2};
  std::vector<double> ts_decreasing = {0.45, 0.2};
  std::vector<double> tsinf = {stan::math::INFTY, 1.1};
  std::vector<double> tsNaN = {0.45, stan::math::NOT_A_NUMBER};

  double a = 1.5;

  std::vector<Eigen::VectorXd> out;
  EXPECT_NO_THROW(out
                  = stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                              1e-10, 1e-10, 1e6, nullptr, a));
  EXPECT_EQ(out.size(), ts.size());

  EXPECT_NO_THROW(out = stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0,
                                                  ts_repeat, 1e-10, 1e-10, 1e6,
                                                  nullptr, a));
  EXPECT_EQ(out.size(), ts_repeat.size());
  EXPECT_MATRIX_FLOAT_EQ(out[0], out[1]);

  EXPECT_NO_THROW(out = stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0,
                                                  ts_lots, 1e-10, 1e-10, 1e6,
                                                  nullptr, a));
  EXPECT_EQ(out.size(), ts_lots.size());

  EXPECT_THROW(
      stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts_empty, 1e-10,
                                1e-10, 1e6, nullptr, a),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts_early, 1e-10,
                                1e-10, 1e6, nullptr, a),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts_decreasing,
                                1e-10, 1e-10, 1e6, nullptr, a),
      std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, tsinf,
                                         1e-10, 1e-10, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, tsNaN,
                                         1e-10, 1e-10, 1e6, nullptr, a),
               std::domain_error);
}

TEST(ode_adams_tol_prim, one_arg_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;
  double ainf = stan::math::INFTY;
  double aNaN = stan::math::NOT_A_NUMBER;

  std::vector<double> va = {a};
  std::vector<double> vainf = {ainf};
  std::vector<double> vaNaN = {aNaN};

  Eigen::VectorXd ea(1);
  ea << a;
  Eigen::VectorXd eainf(1);
  eainf << ainf;
  Eigen::VectorXd eaNaN(1);
  eaNaN << aNaN;

  std::vector<std::vector<double>> vva = {va};
  std::vector<std::vector<double>> vvainf = {vainf};
  std::vector<std::vector<double>> vvaNaN = {vaNaN};

  std::vector<Eigen::VectorXd> vea = {ea};
  std::vector<Eigen::VectorXd> veainf = {eainf};
  std::vector<Eigen::VectorXd> veaNaN = {eaNaN};

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, ainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, aNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, va));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, vainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, vaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, ea));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, eainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, eaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, vva));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, vvainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, vvaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, vea));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, veainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, veaNaN),
               std::domain_error);
}

TEST(ode_adams_tol_prim, two_arg_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;
  double ainf = stan::math::INFTY;
  double aNaN = stan::math::NOT_A_NUMBER;

  std::vector<double> va = {a};
  std::vector<double> vainf = {ainf};
  std::vector<double> vaNaN = {aNaN};

  Eigen::VectorXd ea(1);
  ea << a;
  Eigen::VectorXd eainf(1);
  eainf << ainf;
  Eigen::VectorXd eaNaN(1);
  eaNaN << aNaN;

  std::vector<std::vector<double>> vva = {va};
  std::vector<std::vector<double>> vvainf = {vainf};
  std::vector<std::vector<double>> vvaNaN = {vaNaN};

  std::vector<Eigen::VectorXd> vea = {ea};
  std::vector<Eigen::VectorXd> veainf = {eainf};
  std::vector<Eigen::VectorXd> veaNaN = {eaNaN};

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, a, a));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, ainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, aNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, a, va));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, vainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, vaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                            1e-10, 1e-10, 1e6, nullptr, a, ea));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, eainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, eaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(
      stan::test::Cos2Arg(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a, vva));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, vvainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, vvaNaN),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::ode_adams_tol(
      stan::test::Cos2Arg(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a, vea));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, veainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::Cos2Arg(), y0, t0, ts,
                                         1e-10, 1e-10, 1e6, nullptr, a, veaNaN),
               std::domain_error);
}

TEST(ode_adams_tol_prim, rtol_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double rtol = 1e-6;
  double rtol_negative = -1e-6;
  double rtolinf = stan::math::INFTY;
  double rtolNaN = stan::math::NOT_A_NUMBER;

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            rtol, 1e-10, 1e6, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         rtol_negative, 1e-10, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         rtolinf, 1e-10, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         rtolNaN, 1e-10, 1e6, nullptr, a),
               std::domain_error);
}

TEST(ode_adams_tol_prim, atol_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double atol = 1e-6;
  double atol_negative = -1e-6;
  double atolinf = stan::math::INFTY;
  double atolNaN = stan::math::NOT_A_NUMBER;

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-6, atol, 1e6, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-6, atol_negative, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-6, atolinf, 1e6, nullptr, a),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                         1e-6, atolNaN, 1e6, nullptr, a),
               std::domain_error);
}

TEST(ode_adams_tol_prim, max_num_steps_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  int max_num_steps = 500;
  int max_num_steps_negative = -500;
  int max_num_steps_zero = 0;

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-6, 1e-6, max_num_steps, nullptr,
                                            a));

  EXPECT_THROW(
      stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts, 1e-6, 1e-6,
                                max_num_steps_negative, nullptr, a),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts, 1e-6, 1e-6,
                                max_num_steps_zero, nullptr, a),
      std::domain_error);
}

TEST(ode_adams_tol_prim, rhs_wrong_size_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                            1e-6, 1e-6, 100, nullptr, a));

  EXPECT_THROW(stan::math::ode_adams_tol(stan::test::CosArgWrongSize(), y0, t0,
                                         ts, 1e-6, 1e-6, 100, nullptr, a),
               std::invalid_argument);
}

TEST(ode_adams_tol_prim, error_name) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double ainf = stan::math::INFTY;

  EXPECT_THROW_MSG(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                             1e-6, 1e-6, 100, nullptr, ainf),
                   std::domain_error, "ode_adams_tol");
}

TEST(ode_adams_tol_prim, too_much_work) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1e10};

  double a = 1.0;

  EXPECT_THROW_MSG(stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts,
                                             1e-6, 1e-6, 100, nullptr, a),
                   std::domain_error,
                   "ode_adams_tol:  Failed to integrate to next output time");
}
