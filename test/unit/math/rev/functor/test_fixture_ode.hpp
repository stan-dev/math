#ifndef TEST_UNIT_MATH_REV_FUNCTOR_ODE_FIXTURE_HPP
#define TEST_UNIT_MATH_REV_FUNCTOR_ODE_FIXTURE_HPP

#include <stan/math/rev.hpp>
#include <test/prob/utility.hpp>
#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <type_traits>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

template <class ode_problem_type>
struct ODETestFixture : public ::testing::Test {
  void test_good() {
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);
    ASSERT_NO_THROW(ode.apply_solver());
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

    std::vector<Eigen::VectorXd> results(ode.ts.size());

    for (size_t i = 0; i < ode.ts.size(); ++i) {
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

    std::vector<Eigen::VectorXd> results(ode.ts.size());

    for (size_t i = 0; i < ode.ts.size(); ++i) {
      results[i] = (res_ub[i] - res_lb[i]) / (2.0 * h);
    }
    return results;
  }

  /**
   * Test AD against finite diff when param
   * is <code>var</code>.
   *
   * Require <code>apply_solver(T1&& init, T2&& param)</code> from child fixture
   * for finite diff grad calculation. The call should return ODE data results
   * with when <code>T1</code> and <code>T2</code> are data.
   *
   * @param diff finite diff stepsize
   * @param tol double value test tolerance
   */
  void test_fd_dv(double diff, double tol) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    auto theta_v = stan::math::to_var(ode.param());
    int n = theta_v.size();
    std::vector<std::vector<Eigen::VectorXd>> fd_res(n);
    for (size_t i = 0; i < n; ++i) {
      fd_res[i] = fd_param(i, diff);
    }
    std::vector<double> grads_eff;

    std::vector<Eigen::Matrix<stan::math::var, -1, 1>> ode_res
        = ode.apply_solver(ode.init(), theta_v);

    for (size_t i = 0; i < ode.ts.size(); i++) {
      for (size_t j = 0; j < ode_res[0].size(); j++) {
        grads_eff.clear();
        ode_res[i][j].grad(theta_v, grads_eff);

        for (size_t k = 0; k < n; k++)
          EXPECT_NEAR(grads_eff[k], fd_res[k][i][j], tol)
              << "Gradient of ODE solver failed with initial positions"
              << " known and parameters unknown at time index " << i
              << ", equation index " << j << ", and parameter index: " << k;

        stan::math::set_zero_all_adjoints();
      }
    }
  }

  /**
   * Test AD against finite diff when initial condition
   * is <code>var</code>.
   *
   * Require <code>apply_solver(T1&& init, T2&& param)</code> from child fixture
   * for finite diff grad calculation. The call should return ODE data results
   * with when <code>T1</code> and <code>T2</code> are data.
   *
   * @param diff finite diff stepsize
   * @param tol double value test tolerance
   */
  void test_fd_vd(double diff, double tol) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

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

    for (size_t i = 0; i < ode.ts.size(); i++) {
      for (size_t j = 0; j < n; j++) {
        grads_eff.clear();
        ode_res[i][j].grad(y_vec, grads_eff);

        for (size_t k = 0; k < n; k++)
          EXPECT_NEAR(grads_eff[k], fd_res[k][i][j], tol)
              << "Gradient of ode solver failed with initial positions"
              << " unknown and parameters known at time index " << i
              << ", equation index " << j << ", and parameter index: " << k;

        stan::math::set_zero_all_adjoints();
      }
    }
  }

  /**
   * Test AD against finite diff when both initial condition & param
   * are <code>var</code>.
   *
   * Require <code>apply_solver(T1&& init, T2&& param)</code> from child fixture
   * for finite diff grad calculation. The call should return ODE data results
   * with when <code>T1</code> and <code>T2</code> are data.
   *
   * @param diff finite diff stepsize
   * @param tol double value test tolerance
   */
  void test_fd_vv(double diff, double tol) {
    std::stringstream msgs;
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);

    int n = ode.dim();
    int m = ode.param_size();
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
    for (size_t i = 0; i < ode.ts.size(); i++) {
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

        stan::math::set_zero_all_adjoints();
      }
    }
  }

  /**
   * Test AD when <code>ts</code> is <code>var</code>.
   *
   * require <code>apply_solver()</code> from child fixture for ODE
   * solution when time step is <code>var</code>, and <code>eval_rhs</code>
   * to calculate RHS.
   */
  void test_ts_ad() {
    ode_problem_type& ode = static_cast<ode_problem_type&>(*this);
    std::vector<double> g;
    std::vector<Eigen::Matrix<stan::math::var, -1, 1>> res = ode.apply_solver();
    size_t nt = res.size();
    for (auto i = 0; i < nt; ++i) {
      Eigen::VectorXd res_d = stan::math::value_of(res[i]);
      for (auto j = 0; j < ode.dim(); ++j) {
        g.clear();
        res[i][j].grad();
        for (auto k = 0; k < nt; ++k) {
          if (k != i) {
            EXPECT_FLOAT_EQ(ode.ts[k].adj(), 0.0);
          } else {
            double ts_ad
                = stan::math::value_of(ode.eval_rhs(ode.ts[i].val(), res_d)[j]);
            EXPECT_FLOAT_EQ(ode.ts[k].adj(), ts_ad);
          }
        }
        stan::math::set_zero_all_adjoints();
      }
    }
  }
};

#endif
