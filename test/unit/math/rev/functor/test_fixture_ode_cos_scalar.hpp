#ifndef STAN_MATH_TEST_FIXTURE_ODE_COS_SCALAR_HPP
#define STAN_MATH_TEST_FIXTURE_ODE_COS_SCALAR_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/prim/functor/ode_test_functors.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

struct cos_arg_ode_base {
  stan::test::CosArg1 f1;
  stan::test::Cos2Arg f2;

  Eigen::VectorXd y0;
  double t0;
  std::vector<double> ts;
  double a;
  double rtol;
  double atol;
  int max_num_step;

  cos_arg_ode_base()
      : y0(1),
        t0(0.0),
        ts{0.45, 1.1},
        a(1.5),
        rtol(1.e-10),
        atol(1.e-10),
        max_num_step(100000) {
    y0[0] = 0.0;
  }
};

/**
 * Inheriting base type, various fixtures differs by the type of ODE
 * functor used in <code>apply_solver</code> calls, intended for
 * different kind of tests.
 *
 */
template <typename T>
struct cos_arg_test : public cos_arg_ode_base,
                      public ODETestFixture<cos_arg_test<T>> {
  cos_arg_test() : cos_arg_ode_base() {}

  Eigen::VectorXd init() { return y0; }
  std::vector<double> param() { return {a}; }

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(stan::test::CosArg1(), y0, t0, ts, nullptr, a);
  }

  template <typename T1, typename T2>
  auto apply_solver(Eigen::Matrix<T1, -1, 1>& init, std::vector<T2>& va) {
    std::tuple_element_t<0, T> sol;
    return sol(stan::test::CosArg1(), init, t0, ts, nullptr, va);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(stan::test::CosArg1(), y0, t0, ts, rtol, atol, max_num_step,
               nullptr, a);
  }

  template <typename a_type>
  auto apply_solver_arg(a_type const& a_) {
    std::tuple_element_t<0, T> sol;
    return sol(stan::test::CosArg1(), y0, t0, ts, nullptr, a_);
  }

  template <typename a_type>
  auto apply_solver_arg_tol(a_type const& a_) {
    std::tuple_element_t<1, T> sol;
    return sol(stan::test::CosArg1(), y0, t0, ts, rtol, atol, max_num_step,
               nullptr, a_);
  }

  template <typename a_type, typename b_type>
  auto apply_solver_arg(a_type const& a_, b_type const& b_) {
    std::tuple_element_t<0, T> sol;
    return sol(stan::test::CosArg1(), y0, t0, ts, nullptr, a_, b_);
  }

  template <typename a_type, typename b_type>
  auto apply_solver_arg_tol(a_type const& a_, b_type const& b_) {
    std::tuple_element_t<1, T> sol;
    return sol(stan::test::CosArg1(), y0, t0, ts, rtol, atol, max_num_step,
               nullptr, a_, b_);
  }

  void test_y0_error() {
    y0 = Eigen::VectorXd::Zero(1);
    ASSERT_NO_THROW(apply_solver());

    y0[0] = stan::math::INFTY;
    EXPECT_THROW(apply_solver(), std::domain_error);

    y0[0] = stan::math::NOT_A_NUMBER;
    EXPECT_THROW(apply_solver(), std::domain_error);

    y0 = Eigen::VectorXd();
    EXPECT_THROW(apply_solver(), std::invalid_argument);
  }

  void test_y0_error_with_tol() {
    y0 = Eigen::VectorXd::Zero(1);
    ASSERT_NO_THROW(apply_solver_tol());

    y0[0] = stan::math::INFTY;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    y0[0] = stan::math::NOT_A_NUMBER;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    y0 = Eigen::VectorXd();
    EXPECT_THROW(apply_solver_tol(), std::invalid_argument);
  }

  void test_t0_error() {
    t0 = stan::math::INFTY;
    EXPECT_THROW(apply_solver(), std::domain_error);

    t0 = stan::math::NOT_A_NUMBER;
    EXPECT_THROW(apply_solver(), std::domain_error);
  }

  void test_t0_error_with_tol() {
    t0 = stan::math::INFTY;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    t0 = stan::math::NOT_A_NUMBER;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);
  }

  void test_ts_error() {
    std::vector<double> ts_repeat = {0.45, 0.45};
    std::vector<double> ts_lots = {0.45, 0.45, 1.1, 1.1, 2.0};
    std::vector<double> ts_empty = {};
    std::vector<double> ts_early = {-0.45, 0.2};
    std::vector<double> ts_decreasing = {0.45, 0.2};
    std::vector<double> tsinf = {stan::math::INFTY, 1.1};
    std::vector<double> tsNaN = {0.45, stan::math::NOT_A_NUMBER};

    std::vector<Eigen::VectorXd> out;
    EXPECT_NO_THROW(out = apply_solver());
    EXPECT_EQ(out.size(), ts.size());

    ts = ts_repeat;
    EXPECT_NO_THROW(out = apply_solver());
    EXPECT_EQ(out.size(), ts_repeat.size());
    EXPECT_MATRIX_FLOAT_EQ(out[0], out[1]);

    ts = ts_lots;
    EXPECT_NO_THROW(out = apply_solver());
    EXPECT_EQ(out.size(), ts_lots.size());

    ts = ts_empty;
    EXPECT_THROW(apply_solver(), std::invalid_argument);

    ts = ts_early;
    EXPECT_THROW(apply_solver(), std::domain_error);

    ts = ts_decreasing;
    EXPECT_THROW(apply_solver(), std::domain_error);

    ts = tsinf;
    EXPECT_THROW(apply_solver(), std::domain_error);

    ts = tsNaN;
    EXPECT_THROW(apply_solver(), std::domain_error);

    ts = {0.45, 1.1};
  }

  void test_ts_error_with_tol() {
    std::vector<double> ts_repeat = {0.45, 0.45};
    std::vector<double> ts_lots = {0.45, 0.45, 1.1, 1.1, 2.0};
    std::vector<double> ts_empty = {};
    std::vector<double> ts_early = {-0.45, 0.2};
    std::vector<double> ts_decreasing = {0.45, 0.2};
    std::vector<double> tsinf = {stan::math::INFTY, 1.1};
    std::vector<double> tsNaN = {0.45, stan::math::NOT_A_NUMBER};

    std::vector<Eigen::VectorXd> out;
    EXPECT_NO_THROW(out = apply_solver_tol());
    EXPECT_EQ(out.size(), ts.size());

    ts = ts_repeat;
    EXPECT_NO_THROW(out = apply_solver_tol());
    EXPECT_EQ(out.size(), ts_repeat.size());
    EXPECT_MATRIX_FLOAT_EQ(out[0], out[1]);

    ts = ts_lots;
    EXPECT_NO_THROW(out = apply_solver_tol());
    EXPECT_EQ(out.size(), ts_lots.size());

    ts = ts_empty;
    EXPECT_THROW(apply_solver_tol(), std::invalid_argument);

    ts = ts_early;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    ts = ts_decreasing;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    ts = tsinf;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    ts = tsNaN;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    ts = {0.45, 1.1};
  }

  void test_one_arg_error() {
    a = 1.5;
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

    EXPECT_NO_THROW(apply_solver());

    a = ainf;
    EXPECT_THROW(apply_solver(), std::domain_error);

    a = aNaN;
    EXPECT_THROW(apply_solver(), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(va));
    EXPECT_THROW(apply_solver_arg(vainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg(vaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(ea));
    EXPECT_THROW(apply_solver_arg(eainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg(eaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(vva));
    EXPECT_THROW(apply_solver_arg(vvainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg(vvaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(vea));
    EXPECT_THROW(apply_solver_arg(veainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg(veaNaN), std::domain_error);
  }

  void test_one_arg_error_with_tol() {
    a = 1.5;
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

    EXPECT_NO_THROW(apply_solver_tol());

    a = ainf;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    a = aNaN;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(va));
    EXPECT_THROW(apply_solver_arg_tol(vainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg_tol(vaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(ea));
    EXPECT_THROW(apply_solver_arg_tol(eainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg_tol(eaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(vva));
    EXPECT_THROW(apply_solver_arg_tol(vvainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg_tol(vvaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(vea));
    EXPECT_THROW(apply_solver_arg_tol(veainf), std::domain_error);
    EXPECT_THROW(apply_solver_arg_tol(veaNaN), std::domain_error);
  }

  void test_two_arg_error() {
    a = 1.5;
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

    EXPECT_NO_THROW(apply_solver_arg(a, a));

    EXPECT_THROW(apply_solver_arg(a, ainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg(a, aNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(a, va));

    EXPECT_THROW(apply_solver_arg(a, vainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg(a, vaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(a, ea));

    EXPECT_THROW(apply_solver_arg(a, eainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg(a, eaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(a, vva));

    EXPECT_THROW(apply_solver_arg(a, vvainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg(a, vvaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg(a, vea));

    EXPECT_THROW(apply_solver_arg(a, veainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg(a, veaNaN), std::domain_error);
  }

  void test_two_arg_error_with_tol() {
    a = 1.5;
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

    EXPECT_NO_THROW(apply_solver_arg_tol(a, a));

    EXPECT_THROW(apply_solver_arg_tol(a, ainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg_tol(a, aNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(a, va));

    EXPECT_THROW(apply_solver_arg_tol(a, vainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg_tol(a, vaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(a, ea));

    EXPECT_THROW(apply_solver_arg_tol(a, eainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg_tol(a, eaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(a, vva));

    EXPECT_THROW(apply_solver_arg_tol(a, vvainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg_tol(a, vvaNaN), std::domain_error);

    EXPECT_NO_THROW(apply_solver_arg_tol(a, vea));

    EXPECT_THROW(apply_solver_arg_tol(a, veainf), std::domain_error);

    EXPECT_THROW(apply_solver_arg_tol(a, veaNaN), std::domain_error);
  }

  void test_rhs_wrong_size_error() {
    std::tuple_element_t<0, T> sol;
    EXPECT_THROW(sol(stan::test::CosArgWrongSize(), y0, t0, ts, nullptr, a),
                 std::invalid_argument);
  }

  void test_rhs_wrong_size_error_with_tol() {
    std::tuple_element_t<1, T> sol;
    EXPECT_THROW(sol(stan::test::CosArgWrongSize(), y0, t0, ts, rtol, atol,
                     max_num_step, nullptr, a),
                 std::invalid_argument);
  }

  void test_error_name() {
    a = stan::math::INFTY;
    std::tuple_element_t<0, T> sol;
    EXPECT_THROW_MSG(apply_solver(), std::domain_error, sol.functor_name);
  }

  void test_error_name_with_tol() {
    a = stan::math::INFTY;
    std::tuple_element_t<1, T> sol;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, sol.functor_name);
  }

  void test_rtol_error() {
    y0 = Eigen::VectorXd::Zero(1);
    t0 = 0;
    ts = {0.45, 1.1};
    a = 1.5;

    rtol = 1e-6;
    atol = 1e-6;
    double rtol_negative = -1e-6;
    double rtolinf = stan::math::INFTY;
    double rtolNaN = stan::math::NOT_A_NUMBER;

    EXPECT_NO_THROW(apply_solver_tol());

    rtol = rtol_negative;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    rtol = rtolinf;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    rtol = rtolNaN;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);
  }

  void test_atol_error() {
    y0 = Eigen::VectorXd::Zero(1);
    t0 = 0;
    ts = {0.45, 1.1};
    a = 1.5;

    rtol = 1e-6;
    atol = 1e-6;
    double atol_negative = -1e-6;
    double atolinf = stan::math::INFTY;
    double atolNaN = stan::math::NOT_A_NUMBER;

    EXPECT_NO_THROW(apply_solver_tol());

    atol = atol_negative;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    atol = atolinf;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    atol = atolNaN;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);
  }

  void test_max_num_step_error() {
    rtol = 1e-6;
    atol = 1e-6;
    max_num_step = 500;
    int max_num_steps_negative = -500;
    int max_num_steps_zero = 0;

    EXPECT_NO_THROW(apply_solver_tol());

    max_num_step = max_num_steps_negative;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);

    max_num_step = max_num_steps_zero;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);
  }

  void test_too_much_work() {
    ts[1] = 1e4;
    max_num_step = 10;
    EXPECT_THROW(apply_solver_tol(), std::domain_error);
  }

  void test_value() {
    std::vector<Eigen::VectorXd> res = apply_solver();
    EXPECT_NEAR(res[0][0], 0.4165982112, 1e-5);
    EXPECT_NEAR(res[1][0], 0.66457668563, 1e-5);

    std::vector<double> ts_i = {1, 2};
    std::tuple_element_t<0, T> sol;
    res = sol(stan::test::CosArg1(), y0, t0, ts_i, nullptr, a);
    EXPECT_NEAR(res[0][0], 0.6649966577, 1e-5);
    EXPECT_NEAR(res[1][0], 0.09408000537, 1e-5);

    int t0_i = 0;
    res = sol(stan::test::CosArg1(), y0, t0_i, ts, nullptr, a);
    EXPECT_NEAR(res[0][0], 0.4165982112, 1e-5);
    EXPECT_NEAR(res[1][0], 0.66457668563, 1e-5);
  }

  void test_grad_t0() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    stan::math::var t0v = 0.0;
    auto res = sol(stan::test::CosArg1(), y0, t0v, ts, nullptr, a);

    res[0][0].grad();

    EXPECT_NEAR(res[0][0].val(), 0.4165982112, 1e-5);
    EXPECT_NEAR(t0v.adj(), -1.0, 1e-5);

    nested.set_zero_all_adjoints();

    res[1][0].grad();

    EXPECT_NEAR(res[1][0].val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(t0v.adj(), -1.0, 1e-5);
  }

  void test_grad_ts() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    std::vector<var> tsv = {0.45, 1.1};
    auto res = sol(stan::test::CosArg1(), y0, t0, tsv, nullptr, a);

    res[0][0].grad();

    EXPECT_NEAR(res[0][0].val(), 0.4165982112, 1e-5);
    EXPECT_NEAR(tsv[0].adj(), 0.78070695113, 1e-5);
    nested.set_zero_all_adjoints();

    res[1][0].grad();

    EXPECT_NEAR(res[1][0].val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(tsv[1].adj(), -0.0791208888, 1e-5);
  }

  void test_grad_ts_repeat() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    std::vector<var> tsv = {0.45, 0.45, 1.1, 1.1};
    auto output = sol(stan::test::CosArg1(), y0, t0, tsv, nullptr, a);

    EXPECT_EQ(output.size(), tsv.size());

    output[0][0].grad();

    EXPECT_NEAR(output[0][0].val(), 0.4165982112, 1e-5);
    EXPECT_NEAR(tsv[0].adj(), 0.78070695113, 1e-5);
    nested.set_zero_all_adjoints();

    output[1][0].grad();

    EXPECT_NEAR(output[1][0].val(), 0.4165982112, 1e-5);
    EXPECT_NEAR(tsv[1].adj(), 0.78070695113, 1e-5);
    nested.set_zero_all_adjoints();

    output[2][0].grad();

    EXPECT_NEAR(output[2][0].val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(tsv[2].adj(), -0.0791208888, 1e-5);
    nested.set_zero_all_adjoints();

    output[3][0].grad();
    EXPECT_NEAR(output[3][0].val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(tsv[3].adj(), -0.0791208888, 1e-5);
  }

  void test_scalar_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    stan::math::var av = 1.5;

    {
      std::vector<double> ts1{1.1};
      auto output = sol(stan::test::CosArg1(), y0, t0, ts1, nullptr, av)[0][0];

      output.grad();

      EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
      EXPECT_NEAR(av.adj(), -0.50107310888, 1e-5);
      nested.set_zero_all_adjoints();
    }

    {
      auto output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, av);

      output[0](0).grad();

      EXPECT_NEAR(output[0](0).val(), 0.4165982112, 1e-5);
      EXPECT_NEAR(av.adj(), -0.04352005542, 1e-5);
      nested.set_zero_all_adjoints();

      output[1](0).grad();

      EXPECT_NEAR(output[1](0).val(), 0.66457668563, 1e-5);
      EXPECT_NEAR(av.adj(), -0.50107310888, 1e-5);
      nested.set_zero_all_adjoints();
    }
  }

  void test_std_vector_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    std::vector<var> av = {1.5};
    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, av)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(av[0].adj(), -0.50107310888, 1e-5);
  }

  void test_vector_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    Eigen::Matrix<var, Eigen::Dynamic, 1> av(1);
    av << 1.5;

    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, av)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(av(0).adj(), -0.50107310888, 1e-5);
  }

  void test_row_vector_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    Eigen::Matrix<var, 1, -1> av(1);
    av << 1.5;

    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, av)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(av(0).adj(), -0.50107310888, 1e-5);
  }

  void test_matrix_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    Eigen::Matrix<var, -1, -1> av(1, 1);
    av << 1.5;

    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, av)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(av(0).adj(), -0.50107310888, 1e-5);
  }

  void test_scalar_std_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    var a0 = 0.75;
    std::vector<var> a1 = {0.75};

    var output = sol(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a0, a1)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(a0.adj(), -0.50107310888, 1e-5);
    EXPECT_NEAR(a1[0].adj(), -0.50107310888, 1e-5);
  }

  void test_std_vector_std_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    var a0 = 1.5;
    std::vector<var> a1(1, a0);
    std::vector<std::vector<var>> a2(1, a1);

    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, a2)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(a2[0][0].adj(), -0.50107310888, 1e-5);
  }

  void test_std_vector_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    var a0 = 1.5;
    Eigen::Matrix<var, Eigen::Dynamic, 1> a1(1);
    a1 << a0;
    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> a2(1, a1);

    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, a2)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(a2[0](0).adj(), -0.50107310888, 1e-5);
  }

  void test_std_vector_row_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    var a0 = 1.5;
    Eigen::Matrix<var, 1, Eigen::Dynamic> a1(1);
    a1 << a0;
    std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> a2(1, a1);
    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, a2)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(a2[0](0).adj(), -0.50107310888, 1e-5);
  }

  void test_std_vector_matrix_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;
    var a0 = 1.5;
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a1(1, 1);
    a1 << a0;
    std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> a2(1, a1);

    var output = sol(stan::test::CosArg1(), y0, t0, ts, nullptr, a2)[1][0];

    output.grad();

    EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
    EXPECT_NEAR(a2[0](0).adj(), -0.50107310888, 1e-5);
  }

  void test_arg_combos_test() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<0, T> sol;

    var t0 = 0.5;
    var a = 0.2;
    std::vector<var> ts = {1.25};
    Eigen::Matrix<var, Eigen::Dynamic, 1> y0(1);
    y0 << 0.75;

    double t0d = stan::math::value_of(t0);
    double ad = stan::math::value_of(a);
    std::vector<double> tsd = stan::math::value_of(ts);
    Eigen::VectorXd y0d = stan::math::value_of(y0);

    auto check_yT = [&](auto yT) {
      EXPECT_NEAR(stan::math::value_of(yT),
                  y0d(0) * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
                  1e-5);
    };

    auto check_t0 = [&](var t0) {
      EXPECT_NEAR(
          t0.adj(),
          ad * t0d * y0d(0) * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
          1e-5);
    };

    auto check_a = [&](var a) {
      EXPECT_NEAR(a.adj(),
                  -0.5 * (tsd[0] * tsd[0] - t0d * t0d) * y0d(0)
                      * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
                  1e-5);
    };

    auto check_ts = [&](std::vector<var> ts) {
      EXPECT_NEAR(ts[0].adj(),
                  -ad * tsd[0] * y0d(0)
                      * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
                  1e-5);
    };

    auto check_y0 = [&](Eigen::Matrix<var, Eigen::Dynamic, 1> y0) {
      EXPECT_NEAR(y0(0).adj(), exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
                  1e-5);
    };

    double yT1 = sol(stan::test::ayt(), y0d, t0d, tsd, nullptr, ad)[0](0);
    check_yT(yT1);

    var yT2 = sol(stan::test::ayt(), y0d, t0d, tsd, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT2.grad();
    check_yT(yT2);
    check_a(a);

    var yT3 = sol(stan::test::ayt(), y0d, t0d, ts, nullptr, ad)[0](0);
    nested.set_zero_all_adjoints();
    yT3.grad();
    check_yT(yT3);
    check_ts(ts);

    var yT4 = sol(stan::test::ayt(), y0d, t0d, ts, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT4.grad();
    check_yT(yT4);
    check_ts(ts);
    check_a(a);

    var yT5 = sol(stan::test::ayt(), y0d, t0, tsd, nullptr, ad)[0](0);
    nested.set_zero_all_adjoints();
    yT5.grad();
    check_yT(yT5);
    check_t0(t0);

    var yT6 = sol(stan::test::ayt(), y0d, t0, tsd, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT6.grad();
    check_yT(yT6);
    check_t0(t0);
    check_a(a);

    var yT7 = sol(stan::test::ayt(), y0d, t0, ts, nullptr, ad)[0](0);
    nested.set_zero_all_adjoints();
    yT7.grad();
    check_yT(yT7);
    check_t0(t0);
    check_ts(ts);

    var yT8 = sol(stan::test::ayt(), y0d, t0, ts, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT8.grad();
    check_yT(yT8);
    check_t0(t0);
    check_ts(ts);
    check_a(a);

    var yT9 = sol(stan::test::ayt(), y0, t0d, tsd, nullptr, ad)[0](0);
    nested.set_zero_all_adjoints();
    yT9.grad();
    check_yT(yT9);
    check_y0(y0);

    var yT10 = sol(stan::test::ayt(), y0, t0d, tsd, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT10.grad();
    check_yT(yT10);
    check_y0(y0);
    check_a(a);

    var yT11 = sol(stan::test::ayt(), y0, t0d, ts, nullptr, ad)[0](0);
    nested.set_zero_all_adjoints();
    yT11.grad();
    check_yT(yT11);
    check_y0(y0);
    check_ts(ts);

    var yT12 = sol(stan::test::ayt(), y0, t0d, ts, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT12.grad();
    check_yT(yT12);
    check_y0(y0);
    check_ts(ts);
    check_a(a);

    var yT13 = sol(stan::test::ayt(), y0, t0, tsd, nullptr, ad)[0](0);
    nested.set_zero_all_adjoints();
    yT13.grad();
    check_yT(yT13);
    check_y0(y0);
    check_t0(t0);

    var yT14 = sol(stan::test::ayt(), y0, t0, tsd, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT14.grad();
    check_yT(yT14);
    check_y0(y0);
    check_t0(t0);
    check_a(a);

    var yT15 = sol(stan::test::ayt(), y0, t0, ts, nullptr, ad)[0](0);
    nested.set_zero_all_adjoints();
    yT15.grad();
    check_yT(yT15);
    check_y0(y0);
    check_t0(t0);
    check_ts(ts);

    var yT16 = sol(stan::test::ayt(), y0, t0, ts, nullptr, a)[0](0);
    nested.set_zero_all_adjoints();
    yT16.grad();
    check_yT(yT16);
    check_y0(y0);
    check_t0(t0);
    check_ts(ts);
    check_a(a);
  }

  void test_tol_int_t0() {
    using stan::math::var;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    int t0 = 0;
    std::vector<double> ts = {0.45, 1.1};

    double a = 1.5;

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> output
        = stan::math::ode_adams_tol(stan::test::CosArg1(), y0, t0, ts, 1e-10,
                                    1e-10, 1e6, nullptr, a);

    EXPECT_FLOAT_EQ(output[0][0], 0.4165982112);
    EXPECT_FLOAT_EQ(output[1][0], 0.66457668563);
  }

  void test_tol_int_ts() {
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1, 2};

    double a = 1.5;

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> output
        = sol(stan::test::CosArg1(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a);

    EXPECT_FLOAT_EQ(output[0][0], 0.6649966577);
    EXPECT_FLOAT_EQ(output[1][0], 0.09408000537);
  }

  void test_tol_t0() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    var t0 = 0.0;
    std::vector<double> ts = {0.45, 1.1};

    double a = 1.5;

    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
        = sol(stan::test::CosArg1(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a);

    output[0][0].grad();

    EXPECT_FLOAT_EQ(output[0][0].val(), 0.4165982112);
    EXPECT_FLOAT_EQ(t0.adj(), -1.0);

    nested.set_zero_all_adjoints();

    output[1][0].grad();

    EXPECT_FLOAT_EQ(output[1][0].val(), 0.66457668563);
    EXPECT_FLOAT_EQ(t0.adj(), -1.0);
  }

  void test_tol_ts() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<var> ts = {0.45, 1.1};

    double a = 1.5;

    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
        = sol(stan::test::CosArg1(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a);

    output[0][0].grad();

    EXPECT_FLOAT_EQ(output[0][0].val(), 0.4165982112);
    EXPECT_FLOAT_EQ(ts[0].adj(), 0.78070695113);

    nested.set_zero_all_adjoints();

    output[1][0].grad();

    EXPECT_FLOAT_EQ(output[1][0].val(), 0.66457668563);
    EXPECT_FLOAT_EQ(ts[1].adj(), -0.0791208888);
  }

  void test_tol_ts_repeat() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<var> ts = {0.45, 0.45, 1.1, 1.1};

    double a = 1.5;

    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
        = sol(stan::test::CosArg1(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a);

    EXPECT_EQ(output.size(), ts.size());

    output[0][0].grad();

    EXPECT_FLOAT_EQ(output[0][0].val(), 0.4165982112);
    EXPECT_FLOAT_EQ(ts[0].adj(), 0.78070695113);

    nested.set_zero_all_adjoints();

    output[1][0].grad();

    EXPECT_FLOAT_EQ(output[1][0].val(), 0.4165982112);
    EXPECT_FLOAT_EQ(ts[1].adj(), 0.78070695113);

    nested.set_zero_all_adjoints();

    output[2][0].grad();

    EXPECT_FLOAT_EQ(output[2][0].val(), 0.66457668563);
    EXPECT_FLOAT_EQ(ts[2].adj(), -0.0791208888);

    nested.set_zero_all_adjoints();

    output[3][0].grad();

    EXPECT_FLOAT_EQ(output[3][0].val(), 0.66457668563);
    EXPECT_FLOAT_EQ(ts[3].adj(), -0.0791208888);
  }

  void test_tol_scalar_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    var a = 1.5;

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a.adj(), -0.50107310888);
  }

  void test_tol_scalar_arg_multi_time() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {0.45, 1.1};

    var a = 1.5;

    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
        = sol(stan::test::CosArg1(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a);

    output[0](0).grad();

    EXPECT_FLOAT_EQ(output[0](0).val(), 0.4165982112);
    EXPECT_FLOAT_EQ(a.adj(), -0.04352005542);

    nested.set_zero_all_adjoints();

    output[1](0).grad();

    EXPECT_FLOAT_EQ(output[1](0).val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a.adj(), -0.50107310888);
  }

  void test_tol_std_vector_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    std::vector<var> a = {1.5};

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a[0].adj(), -0.50107310888);
  }

  void test_tol_vector_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    Eigen::Matrix<var, Eigen::Dynamic, 1> a(1);
    a << 1.5;

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a(0).adj(), -0.50107310888);
  }

  void test_tol_row_vector_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    Eigen::Matrix<var, 1, Eigen::Dynamic> a(1);
    a << 1.5;

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a(0).adj(), -0.50107310888);
  }

  void test_tol_matrix_arg() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a(1, 1);
    a << 1.5;

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a(0, 0).adj(), -0.50107310888);
  }

  void test_tol_scalar_std_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    var a0 = 0.75;
    std::vector<var> a1 = {0.75};

    var output = sol(stan::test::Cos2Arg(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a0, a1)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a0.adj(), -0.50107310888);
    EXPECT_FLOAT_EQ(a1[0].adj(), -0.50107310888);
  }

  void test_tol_std_vector_std_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    var a0 = 1.5;
    std::vector<var> a1(1, a0);
    std::vector<std::vector<var>> a2(1, a1);

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a2)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a2[0][0].adj(), -0.50107310888);
  }

  void test_tol_std_vector_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    var a0 = 1.5;
    Eigen::Matrix<var, Eigen::Dynamic, 1> a1(1);
    a1 << a0;
    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> a2(1, a1);

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a2)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a2[0](0).adj(), -0.50107310888);
  }

  void test_tol_std_vector_row_vector_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    var a0 = 1.5;
    Eigen::Matrix<var, 1, Eigen::Dynamic> a1(1);
    a1 << a0;
    std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> a2(1, a1);

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a2)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a2[0](0).adj(), -0.50107310888);
  }

  void test_tol_std_vector_matrix_args() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
    double t0 = 0.0;
    std::vector<double> ts = {1.1};

    var a0 = 1.5;
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a1(1, 1);
    a1 << a0;
    std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> a2(1, a1);

    var output = sol(stan::test::CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                     nullptr, a2)[0][0];

    output.grad();

    EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
    EXPECT_FLOAT_EQ(a2[0](0).adj(), -0.50107310888);
  }

  void test_tol_arg_combos_test() {
    stan::math::nested_rev_autodiff nested;
    std::tuple_element_t<1, T> sol;

    var t0 = 0.5;
    var a = 0.2;
    std::vector<var> ts = {1.25};
    Eigen::Matrix<var, Eigen::Dynamic, 1> y0(1);
    y0 << 0.75;

    double t0d = stan::math::value_of(t0);
    double ad = stan::math::value_of(a);
    std::vector<double> tsd = stan::math::value_of(ts);
    Eigen::VectorXd y0d = stan::math::value_of(y0);

    auto check_yT = [&](auto yT) {
      EXPECT_FLOAT_EQ(stan::math::value_of(yT),
                      y0d(0) * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)));
    };

    auto check_t0 = [&](var t0) {
      EXPECT_FLOAT_EQ(
          t0.adj(),
          ad * t0d * y0d(0) * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)));
    };

    auto check_a = [&](var a) {
      EXPECT_FLOAT_EQ(a.adj(),
                      -0.5 * (tsd[0] * tsd[0] - t0d * t0d) * y0d(0)
                          * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)));
    };

    auto check_ts = [&](std::vector<var> ts) {
      EXPECT_FLOAT_EQ(ts[0].adj(),
                      -ad * tsd[0] * y0d(0)
                          * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)));
    };

    auto check_y0 = [&](Eigen::Matrix<var, Eigen::Dynamic, 1> y0) {
      EXPECT_FLOAT_EQ(y0(0).adj(),
                      exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)));
    };

    double yT1 = sol(stan::test::ayt(), y0d, t0d, tsd, 1e-10, 1e-10, 1e6,
                     nullptr, ad)[0](0);
    check_yT(yT1);

    var yT2 = sol(stan::test::ayt(), y0d, t0d, tsd, 1e-10, 1e-10, 1e6, nullptr,
                  a)[0](0);
    nested.set_zero_all_adjoints();
    yT2.grad();
    check_yT(yT2);
    check_a(a);

    var yT3 = sol(stan::test::ayt(), y0d, t0d, ts, 1e-10, 1e-10, 1e6, nullptr,
                  ad)[0](0);
    nested.set_zero_all_adjoints();
    yT3.grad();
    check_yT(yT3);
    check_ts(ts);

    var yT4 = sol(stan::test::ayt(), y0d, t0d, ts, 1e-10, 1e-10, 1e6, nullptr,
                  a)[0](0);
    nested.set_zero_all_adjoints();
    yT4.grad();
    check_yT(yT4);
    check_ts(ts);
    check_a(a);

    var yT5 = sol(stan::test::ayt(), y0d, t0, tsd, 1e-10, 1e-10, 1e6, nullptr,
                  ad)[0](0);
    nested.set_zero_all_adjoints();
    yT5.grad();
    check_yT(yT5);
    check_t0(t0);

    var yT6 = sol(stan::test::ayt(), y0d, t0, tsd, 1e-10, 1e-10, 1e6, nullptr,
                  a)[0](0);
    nested.set_zero_all_adjoints();
    yT6.grad();
    check_yT(yT6);
    check_t0(t0);
    check_a(a);

    var yT7 = sol(stan::test::ayt(), y0d, t0, ts, 1e-10, 1e-10, 1e6, nullptr,
                  ad)[0](0);
    nested.set_zero_all_adjoints();
    yT7.grad();
    check_yT(yT7);
    check_t0(t0);
    check_ts(ts);

    var yT8 = sol(stan::test::ayt(), y0d, t0, ts, 1e-10, 1e-10, 1e6, nullptr,
                  a)[0](0);
    nested.set_zero_all_adjoints();
    yT8.grad();
    check_yT(yT8);
    check_t0(t0);
    check_ts(ts);
    check_a(a);

    var yT9 = sol(stan::test::ayt(), y0, t0d, tsd, 1e-10, 1e-10, 1e6, nullptr,
                  ad)[0](0);
    nested.set_zero_all_adjoints();
    yT9.grad();
    check_yT(yT9);
    check_y0(y0);

    var yT10 = sol(stan::test::ayt(), y0, t0d, tsd, 1e-10, 1e-10, 1e6, nullptr,
                   a)[0](0);
    nested.set_zero_all_adjoints();
    yT10.grad();
    check_yT(yT10);
    check_y0(y0);
    check_a(a);

    var yT11 = sol(stan::test::ayt(), y0, t0d, ts, 1e-10, 1e-10, 1e6, nullptr,
                   ad)[0](0);
    nested.set_zero_all_adjoints();
    yT11.grad();
    check_yT(yT11);
    check_y0(y0);
    check_ts(ts);

    var yT12 = sol(stan::test::ayt(), y0, t0d, ts, 1e-10, 1e-10, 1e6, nullptr,
                   a)[0](0);
    nested.set_zero_all_adjoints();
    yT12.grad();
    check_yT(yT12);
    check_y0(y0);
    check_ts(ts);
    check_a(a);

    var yT13 = sol(stan::test::ayt(), y0, t0, tsd, 1e-10, 1e-10, 1e6, nullptr,
                   ad)[0](0);
    nested.set_zero_all_adjoints();
    yT13.grad();
    check_yT(yT13);
    check_y0(y0);
    check_t0(t0);

    var yT14 = sol(stan::test::ayt(), y0, t0, tsd, 1e-10, 1e-10, 1e6, nullptr,
                   a)[0](0);
    nested.set_zero_all_adjoints();
    yT14.grad();
    check_yT(yT14);
    check_y0(y0);
    check_t0(t0);
    check_a(a);

    var yT15 = sol(stan::test::ayt(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr,
                   ad)[0](0);
    nested.set_zero_all_adjoints();
    yT15.grad();
    check_yT(yT15);
    check_y0(y0);
    check_t0(t0);
    check_ts(ts);

    var yT16 = sol(stan::test::ayt(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr,
                   a)[0](0);
    nested.set_zero_all_adjoints();
    yT16.grad();
    check_yT(yT16);
    check_y0(y0);
    check_t0(t0);
    check_ts(ts);
    check_a(a);
  }
};

#endif
