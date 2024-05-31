#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/solve_newton.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/functor/util_algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//////////////////////////////////////////////////////////////////////////
// Tests for newton solver.

TEST_F(algebra_solver_simple_eq_test, newton_dbl) {
  bool is_newton = true;
  Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y_dbl, is_newton);
}

TEST_F(algebra_solver_simple_eq_test, newton_tuned_dbl) {
  bool is_newton = true;
  Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y_dbl, is_newton,
                                         true, scale_step, xtol, ftol, maxfev);
}

TEST_F(algebra_solver_simple_eq_nopara_test, newton_dbl) {
  using stan::math::solve_newton;
  Eigen::VectorXd theta = solve_newton(simple_eq_functor_nopara(), x,
                                       &std::cout, y_dummy, dat, dummy_dat_int);
  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
}

TEST_F(algebra_solver_non_linear_eq_test, newton_dbl) {
  bool is_newton = true;
  Eigen::VectorXd theta
      = non_linear_eq_test(non_linear_eq_functor(), y_dbl, is_newton);
  EXPECT_FLOAT_EQ(-y_dbl(0), theta(0));
  EXPECT_FLOAT_EQ(-y_dbl(1), theta(1));
  EXPECT_FLOAT_EQ(y_dbl(2), theta(2));
}

TEST_F(error_message_test, newton_dbl) {
  bool is_newton = true;
  error_conditions_test(non_linear_eq_functor(), y_3, is_newton);
}

TEST_F(max_steps_test, newton_dbl) {
  bool is_newton = true;
  max_num_steps_test(y, is_newton);
}

TEST(MathMatrixRevMat, unsolvable_flag_newton_dbl) {
  Eigen::VectorXd y(2);
  y << 1, 1;

  unsolvable_flag_test(y, 1);
}

TEST_F(degenerate_eq_test, newton_guess1_dbl) {
  using stan::math::solve_newton;

  // This first initial guess produces the
  // solution x = {8, 8}

  Eigen::VectorXd theta = solve_newton(degenerate_eq_functor(), x_guess_1,
                                       &std::cout, y_dbl, dat, dat_int);
  EXPECT_FLOAT_EQ(8, theta(0));
  EXPECT_FLOAT_EQ(8, theta(1));
}

TEST_F(degenerate_eq_test, newton_guess2_dbl) {
  using stan::math::solve_newton;
  // This next initial guess produces the
  // solution x = {5, 5}

  Eigen::VectorXd theta = solve_newton(degenerate_eq_functor(), x_guess_2,
                                       &std::cout, y_dbl, dat, dat_int);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));
}

// For the next two unit tests,  see if the initial
// guess determines neighborhood of the
// solution, when solutions have different scales,
// using y_scale.

TEST_F(degenerate_eq_test, newton_guess2_scale_dbl) {
  using stan::math::solve_newton;

  Eigen::VectorXd theta = solve_newton(degenerate_eq_functor(), x_guess_2,
                                       &std::cout, y_scale, dat, dat_int);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));
}

TEST_F(degenerate_eq_test, newton_guess_saddle_point_dbl) {
  // Newton solver fails this test because the initial point is
  // a saddle point.
  using stan::math::solve_newton;
  std::stringstream err_msg;
  err_msg << "The linear solver’s setup function failed in an unrecoverable "
             "manner";  // NOLINT
  std::string msg = err_msg.str();

  EXPECT_THROW_MSG(solve_newton(degenerate_eq_functor(), x_guess_3, &std::cout,
                                y_scale, dat, dat_int),
                   std::runtime_error, msg);
}

TEST_F(algebra_solver_simple_eq_test, newton) {
  using stan::math::var;
  bool is_newton = true;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y, is_newton);

    std::vector<stan::math::var> y_vec{y(0), y(1), y(2)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

TEST_F(algebra_solver_simple_eq_test, newton_tuned) {
  using stan::math::var;
  bool is_newton = true;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y, is_newton, true, scale_step,
                         xtol, ftol, maxfev);

    std::vector<stan::math::var> y_vec{y(0), y(1), y(2)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

TEST_F(algebra_solver_simple_eq_test, newton_init_is_para) {
  using stan::math::solve_newton;
  Eigen::VectorXd theta = solve_newton(simple_eq_functor(), x_var, &std::cout,
                                       y_dbl, dat, dat_int);
  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
}

TEST_F(algebra_solver_non_linear_eq_test, newton) {
  using stan::math::var;
  bool is_newton = true;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = non_linear_eq_test(non_linear_eq_functor(), y, is_newton);

    EXPECT_FLOAT_EQ(-y(0).val(), theta(0).val());
    EXPECT_FLOAT_EQ(-y(1).val(), theta(1).val());
    EXPECT_FLOAT_EQ(y(2).val(), theta(2).val());

    std::vector<stan::math::var> y_vec{y(0), y(1), y(2)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_NEAR(J(k, i), g[i], err);
  }
}

TEST_F(error_message_test, newton) {
  using stan::math::var;
  bool is_newton = true;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_2;
  error_conditions_test(non_linear_eq_functor(), y, is_newton);
}

TEST_F(max_steps_test, newton) {
  bool is_newton = true;
  max_num_steps_test(y_var, is_newton);
}

TEST(MathMatrixRevMat, unsolvable_flag_newton) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y(2);
  y << 1, 1;

  unsolvable_flag_test(y, 1);
}

TEST_F(degenerate_eq_test, newton_guess1) {
  using stan::math::solve_newton;
  using stan::math::var;

  // This first initial guess produces the solution x = {8, 8}.
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = solve_newton(
        degenerate_eq_functor(), x_guess_1, &std::cout, y, dat, dat_int);
    EXPECT_FLOAT_EQ(8, theta(0).val());
    EXPECT_FLOAT_EQ(8, theta(1).val());

    std::vector<stan::math::var> y_vec{y(0), y(1)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J1(k, l), g[l], tolerance);
  }
}

TEST_F(degenerate_eq_test, newton_guess2) {
  using stan::math::solve_newton;
  using stan::math::var;
  // This next initial guess produces the solution x = {5, 5}.
  for (int k = 0; k < 1; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = solve_newton(
        degenerate_eq_functor(), x_guess_2, &std::cout, y, dat, dat_int);
    EXPECT_FLOAT_EQ(5, theta(0).val());
    EXPECT_FLOAT_EQ(5, theta(0).val());

    std::vector<stan::math::var> y_vec{y(0), y(1)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J2(k, l), g[l], tolerance);
  }
}

TEST_F(variadic_test, newton) {
  using stan::math::var;
  bool is_newton = true;
  bool use_tol = false;
  for (int k = 0; k < n_x; k++) {
    var y_1 = y_1_dbl;
    var y_2 = y_2_dbl;
    var y_3 = y_3_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = variadic_eq_impl_test(
        A, y_1, y_2, y_3, i, is_newton, use_tol, scaling_step_size,
        relative_tolerance, function_tolerance, max_num_steps);
    std::vector<var> y_vec{y_1, y_2, y_3};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_NEAR(J(k, i), g[i], 1e-6);
  }
}

// Additional tests for deprecated signature (with and without tol)
TEST_F(algebra_solver_simple_eq_test, newton_deprecated) {
  using stan::math::var;
  bool is_newton = true;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_non_varia_test(simple_eq_non_varia_functor(), y, is_newton);

    std::vector<stan::math::var> y_vec{y(0), y(1), y(2)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

TEST_F(algebra_solver_simple_eq_test, newton_tuned_deprecated) {
  using stan::math::var;
  bool is_newton = true;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_non_varia_test(simple_eq_non_varia_functor(), y, is_newton,
                                   true, scale_step, xtol, ftol, maxfev);

    std::vector<stan::math::var> y_vec{y(0), y(1), y(2)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}
