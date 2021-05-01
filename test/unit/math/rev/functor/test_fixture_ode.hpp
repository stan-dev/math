#ifndef TEST_UNIT_MATH_REV_FUNCTOR_ODE_FIXTURE_HPP
#define TEST_UNIT_MATH_REV_FUNCTOR_ODE_FIXTURE_HPP

#include <stan/math/rev.hpp>
#include <test/prob/utility.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <type_traits>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

/**
 * The parent test fixture for ODEs, with CRTP dependence on child
 * fixtures. The goal is to provide basic tests for various integrator
 * setups without losing flexibility in adding dedicated tests in the future.
 * In order to focus on actual tests, future maintenance
 * should avoid overengineer this fixture.
 *
 * In order to test various ODE input types, child fixture can use
 * tuple pattern + googletest's typed tests
 *
 * https://github.com/google/googletest/blob/master/googletest/samples/sample6_unittest.cc
 *
 * -------------------------
 * Write child test fixtures
 * -------------------------
 *
 * template<class tuple_type>
 * struct foo_ode_test : public ODETestFixture<foo_ode_test<T>> {
 *   using ode_solver_type = std::tuple_element_t<0, T>;
 *   using T_1 = std::tuple_element_t<1, T>;
 *   using T_2 = std::tuple_element_t<2, T>;
 *   using T_3 = std::tuple_element_t<3, T>;
 * ...
 * };
 *
 * Functor wrappers for Math's ODE solvers can be found in
 * "ode_test_functors.hpp". Currently we have
 *
 * - ode_adams_functor
 * - ode_ckrk_functor
 * - ode_bdf_functor
 * - ode_rk45_functor
 * - ode_adjoint_functor
 * - integrate_ode_adams_functor
 * - integrate_ode_bdf_functor
 * - integrate_ode_rk45_functor
 *
 * corresponding to existing integrator functions, with same calling
 * signatures that supports w & w/o tolerance controls.
 *
 * This fixture provides the following test methods:
 *
 *  | Method          | What's being tested                                 |
 *  |-----------------+-----------------------------------------------------|
 *  | test_good       | Passes without throwing exceptions                  |
 *  | test_analytical | matches analytical soln                             |
 *  | test_fd_dv      | matches finite difference soln (theta is var)       |
 *  | test_fd_vd      | matches finite difference soln (y0 is var)          |
 *  | test_fd_vv      | matches finite difference soln (y0 & theta are var) |
 *  | test_ts_ad      | matches finite difference soln (time is var)        |
 *
 * See each method's description for functions required in the child
 * fixture, and note that child fixture only need to provide required
 * methods should it use a particular test.
 *
 * Since ODE tests require specific compile instructions, the test
 * files should be consistently named as follow:
 *
 * | child fixture            | tests                         |
 * |--------------------------+-------------------------------|
 * | test_fixture_ode_foo.hpp | foo_ode_typed_bar_test.cpp    |
 * | e.g.                     | e.g.                          |
 * | test_fixture_ode_sho.hpp | sho_ode_typed_test.cpp,       |
 * |                          | sho_ode_typed_error_test.cpp, |
 * |                          | sho_ode_typed_...             |
 *
 *
 * ----------------------------
 * Use fixtures in tests
 * ----------------------------
 * Follow googletest manual, first we need to declare the test fixture
 *
 *  TYPED_TEST_SUITE_P(foo_test);
 *
 * then we add test items
 *
 *  TYPED_TEST_P(foo_test, this) { this->test_this(); }
 *  TYPED_TEST_P(foo_test, that) { this->test_that(); }
 *  TYPED_TEST_P(foo_test, another) { this->test_another(); }
 *
 * and register test items
 *
 *  REGISTER_TYPED_TEST_SUITE_P(foo_test, this, that, another);
 *
 * To run actual tests, we need to declare types to be tested and
 *  apply it to the tests
 *
 * using foot_test_types = ::testing::Types<std::tuple<T1, T2, T3>,
 *                                          std::tuple<T4, T5, T6>,
 *                                          ... >
 *  INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, foo_test, foo_test_types);
 *
 *
 * For a simple exmaple, see googletest link above &
 * "fho_ode_typed_ts_test.cpp".
 */
template <class ode_problem_type>
struct ODETestFixture : public ::testing::Test {
  virtual void TearDown() { stan::math::recover_memory(); }

  /**
   * test ODE solver pass
   *
   * Require method in child fixture:
   * - apply_solver(): call ODE solver
   */
  void test_good() {
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);
    ASSERT_NO_THROW(ode.apply_solver());
  }

  /**
   * test ODE solution against analytical solution
   *
   * Require functors in args:
   * - Matrix<T, -1, 1> F_ode(VectorXd): return ODE solution given
   *   parameter vector
   * - Matrix<double, -1, 1> F_sol(double t, ...): return analytical solution at
   * t
   *
   * @param ode_sol solver functor that takes a <code>vector</code>
   * parameter variable that returns solution <code>vector</code>.
   * @param analy_sol analytical solution functor that returns
   * solution <code>vector</code>.
   * @param tol comparison tolerance
   * @param x parameter <code>vector</code> fed into <code>ode_sol</code>
   * @param t time at which analyitical solution is evaluated
   * @param args parameters pack required by <code>analy_sol</code>.
   */
  template <typename F_ode, typename F_sol, typename... T_args>
  void test_analytical(F_ode const& ode_sol, F_sol const& analy_sol, double tol,
                       Eigen::VectorXd const& x, double t,
                       const T_args&... args) {
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);
    auto sol = ode_sol(x);
    int n = sol.size();
    auto sol_0 = analy_sol(t, args...);
    for (int i = 0; i < n; ++i) {
      EXPECT_NEAR(stan::math::value_of(sol[i]), sol_0[i], tol);
    }
  }

  /**
   * test ODE solution as well as sensitivity solution
   * against analytical solution.
   *
   *
   * Require functors in child fixture:
   * - Matrix<T, -1, 1> F_ode(vector): return ODE solution given parameter
   * - Matrix<double, -1, 1> F_sol(double t, ...): return analytical ODE
   * solution at t
   * - Matrix<double, -1, -1> F_grad_sol(double t, ...): return analytical ODE
   * sensivity at t
   *
   * @param ode_sol solver functor that takes a <code>vector</code>
   * parameter variable that returns solution <code>vector</code>.
   * @param analy_sol analytical solution functor that returns
   * solution <code>vector</code>.
   * @param analy_grad_sol analytical sensitivity solution functor
   * that returns sensitivity <code>matrix</code>, with column i
   * corresponding to sensitivity of state i w.r.t parameters
   * <code>vector</code> x.
   * @param tol comparison tolerance
   * @param x parameter <code>vector</code> fed into <code>ode_sol</code>
   * @param t time at which analyitical solution is evaluated
   * @param args parameters pack required by <code>analy_sol</code>.
   */
  template <typename F_ode, typename F_sol, typename F_grad_sol,
            typename... T_args>
  void test_analytical(F_ode const& ode_sol, F_sol const& analy_sol,
                       F_grad_sol const& analy_grad_sol, double tol,
                       Eigen::VectorXd const& x, double t,
                       const T_args&... args) {
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    stan::math::nested_rev_autodiff nested;

    auto sol_0 = analy_sol(t, args...);
    Eigen::Matrix<var, -1, 1> x_var(stan::math::to_var(x));
    Eigen::Matrix<stan::math::var, -1, 1> sol = ode_sol(x_var);
    EXPECT_TRUE(sol.size() == sol_0.size());
    for (auto i = 0; i < sol.size(); ++i) {
      EXPECT_NEAR(sol[i].val(), sol_0[i], tol)
          << "ODE solution failed for state i, i = " << i;
    }

    Eigen::Matrix<double, -1, -1> grad_0(analy_grad_sol(t, args...));
    for (auto i = 0; i < sol.size(); ++i) {
      nested.set_zero_all_adjoints();
      sol[i].grad();
      for (auto j = 0; j < x.size(); ++j) {
        EXPECT_NEAR(x_var[j].adj(), grad_0(j, i), tol)
            << "ODE sensitivity solution failed for state i and parameter j, "
               "(i, j) = ("
            << i << ", " << j << ").";
      }
    }
  }

  /**
   * Gradient wrt to certain param using central difference
   *
   * @param param_index index to param of which sensitivity is seeked
   * @param h finite difference step/perturbation.
   *
   * @return gradient
   */
  std::vector<Eigen::VectorXd> fd_param(const size_t& param_index,
                                        const double& h) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    auto param = ode.param();
    auto init = ode.init();
    param[param_index] += h;
    std::vector<Eigen::VectorXd> res_ub = ode.apply_solver(init, param);
    param[param_index] -= 2 * h;
    std::vector<Eigen::VectorXd> res_lb = ode.apply_solver(init, param);

    std::vector<Eigen::VectorXd> results(ode.times().size());

    for (size_t i = 0; i < ode.times().size(); ++i) {
      results[i] = (res_ub[i] - res_lb[i]) / (2.0 * h);
    }
    return results;
  }

  /**
   * Gradient wrt to certain param using central difference
   *
   * @param param_index index to param of which sensitivity is seeked
   * @param h finite difference step/perturbation.
   *
   * @return gradient
   */
  std::vector<Eigen::VectorXd> fd_init(const size_t& param_index,
                                       const double& h) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    auto param = ode.param();
    auto init = ode.init();
    init[param_index] += h;
    std::vector<Eigen::VectorXd> res_ub = ode.apply_solver(init, param);
    init[param_index] -= 2 * h;
    std::vector<Eigen::VectorXd> res_lb = ode.apply_solver(init, param);

    std::vector<Eigen::VectorXd> results(ode.times().size());

    for (size_t i = 0; i < ode.times().size(); ++i) {
      results[i] = (res_ub[i] - res_lb[i]) / (2.0 * h);
    }
    return results;
  }

  /**
   * Test AD against finite diff when param
   * is <code>var</code>.
   *
   * Require methods in child fixture:
   * - param(): return to ODE vector-like parameters
   * - init(): return ODE vector-like init condition
   * - times(): return ODE vector time step
   * - apply_solver(init, param) solve ode given vector-like init & param.
   *   It must return data results when both init & param are data.
   *
   * @param diff finite diff stepsize
   * @param tol double value test tolerance
   */
  void test_fd_dv(double diff, double tol) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    stan::math::nested_rev_autodiff nested;

    auto theta_v = stan::math::to_var(ode.param());
    int n = theta_v.size();
    std::vector<std::vector<Eigen::VectorXd>> fd_res(n);
    for (size_t i = 0; i < n; ++i) {
      fd_res[i] = fd_param(i, diff);
    }
    std::vector<double> grads_eff;

    std::vector<Eigen::Matrix<stan::math::var, -1, 1>> ode_res
        = ode.apply_solver(ode.init(), theta_v);

    for (size_t i = 0; i < ode.times().size(); i++) {
      for (size_t j = 0; j < ode_res[0].size(); j++) {
        grads_eff.clear();
        ode_res[i][j].grad(theta_v, grads_eff);

        for (size_t k = 0; k < n; k++)
          EXPECT_NEAR(grads_eff[k], fd_res[k][i][j], tol)
              << "Gradient of ODE solver failed with initial positions"
              << " known and parameters unknown at time index " << i
              << ", equation index " << j << ", and parameter index: " << k;

        nested.set_zero_all_adjoints();
      }
    }
  }

  /**
   * Test AD against finite diff when initial condition
   * is <code>var</code>.
   *
   * Require methods in child fixture:
   * - param(): return to ODE vector-like parameters
   * - init(): return ODE vector-like init condition
   * - times(): return ODE vector time step
   * - apply_solver(init, param) solve ode given vector-like init & param.
   *   It must return data results when both init & param are data.
   *
   * @param diff finite diff stepsize
   * @param tol double value test tolerance
   */
  void test_fd_vd(double diff, double tol) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    stan::math::nested_rev_autodiff nested;

    Eigen::Matrix<stan::math::var, -1, 1> y0_v = stan::math::to_var(ode.init());
    int n = y0_v.size();
    std::vector<std::vector<Eigen::VectorXd>> fd_res(n);
    for (size_t i = 0; i < n; ++i) {
      fd_res[i] = fd_init(i, diff);
    }
    std::vector<double> grads_eff;

    std::vector<Eigen::Matrix<stan::math::var, -1, 1>> ode_res
        = ode.apply_solver(y0_v, ode.param());

    std::vector<stan::math::var> y_vec(to_array_1d(y0_v));

    for (size_t i = 0; i < ode.times().size(); i++) {
      for (size_t j = 0; j < n; j++) {
        grads_eff.clear();
        ode_res[i][j].grad(y_vec, grads_eff);

        for (size_t k = 0; k < n; k++)
          EXPECT_NEAR(grads_eff[k], fd_res[k][i][j], tol)
              << "Gradient of ode solver failed with initial positions"
              << " unknown and parameters known at time index " << i
              << ", equation index " << j << ", and parameter index: " << k;

        nested.set_zero_all_adjoints();
      }
    }
  }

  /**
   * Test AD against finite diff when both initial condition & param
   * are <code>var</code>.
   *
   * Require methods in child fixture:
   * - param(): return to ODE vector-like parameters
   * - init(): return ODE vector-like init condition
   * - times(): return ODE vector time step
   * - apply_solver(init, param) solve ode given vector-like init & param.
   *   It must return data results when both init & param are data.
   *
   * @param diff finite diff stepsize
   * @param tol double value test tolerance
   */
  void test_fd_vv(double diff, double tol) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    stan::math::nested_rev_autodiff nested;

    int n = ode.init().size();
    int m = ode.param().size();
    std::vector<std::vector<Eigen::VectorXd>> fd_res_y(n);
    for (size_t i = 0; i < n; ++i) {
      fd_res_y[i] = fd_init(i, diff);
    }

    std::vector<std::vector<Eigen::VectorXd>> fd_res_p(m);
    for (size_t i = 0; i < m; ++i) {
      fd_res_p[i] = fd_param(i, diff);
    }

    std::vector<stan::math::var> vars(
        stan::math::to_array_1d(stan::math::to_var(ode.init())));
    auto theta = ode.param();
    for (int i = 0; i < m; ++i) {
      vars.push_back(theta[i]);
    }
    Eigen::Matrix<stan::math::var, -1, 1> yv(n);
    for (int i = 0; i < n; ++i) {
      yv(i) = vars[i];
    }
    std::vector<stan::math::var> theta_v(vars.begin() + n, vars.end());

    std::vector<Eigen::Matrix<stan::math::var, -1, 1>> ode_res
        = ode.apply_solver(yv, theta_v);

    std::vector<double> grads_eff;
    for (size_t i = 0; i < ode.times().size(); i++) {
      for (size_t j = 0; j < n; j++) {
        grads_eff.clear();
        ode_res[i][j].grad(vars, grads_eff);

        for (size_t k = 0; k < m; k++) {
          EXPECT_NEAR(grads_eff[k + n], fd_res_p[k][i][j], tol)
              << "Gradient of ode solver failed with initial positions"
              << " unknown and parameters unknown for param at time index " << i
              << ", equation index " << j << ", and parameter index: " << k;
        }
        for (size_t k = 0; k < n; k++) {
          EXPECT_NEAR(grads_eff[k], fd_res_y[k][i][j], tol)
              << "Gradient of ode solver failed with initial positions"
              << " unknown and parameters known for initial position at time "
                 "index "
              << i << ", equation index " << j
              << ", and parameter index: " << k;
        }

        nested.set_zero_all_adjoints();
      }
    }
  }

  /**
   * Test gradient w.r.t times.
   *
   * Require methods in child fixture:
   * - times(): return ODE var vector time steps
   * - apply_solver(): solve ODE assuming var time step
   * - eval_rhs(double t, VectorXd x): eval ODE right-hand-side at
   *   time t and independent variable x.
   */
  void test_ts_ad() {
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);
    std::vector<double> g;
    std::vector<Eigen::Matrix<stan::math::var, -1, 1>> res = ode.apply_solver();
    size_t nt = res.size();
    for (auto i = 0; i < nt; ++i) {
      Eigen::VectorXd res_d = stan::math::value_of(res[i]);
      for (auto j = 0; j < ode.init().size(); ++j) {
        g.clear();
        res[i][j].grad();
        for (auto k = 0; k < nt; ++k) {
          if (k != i) {
            EXPECT_FLOAT_EQ(ode.times()[k].adj(), 0.0);
          } else {
            double ts_ad = stan::math::value_of(
                ode.eval_rhs(ode.times()[i].val(), res_d)[j]);
            EXPECT_FLOAT_EQ(ode.times()[k].adj(), ts_ad);
          }
        }
        stan::math::set_zero_all_adjoints();
      }
    }
  }
};

#endif
