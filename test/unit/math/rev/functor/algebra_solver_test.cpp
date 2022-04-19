#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_solver_powell.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/functor/util_algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//////////////////////////////////////////////////////////////////////////
// Tests for powell solver.

TEST_F(algebra_solver_simple_eq_test, powell_dbl) {
  bool is_newton = false;
  Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y_dbl, is_newton);
}

TEST_F(algebra_solver_simple_eq_test, powell_tuned_dbl) {
  bool is_newton = false;
  Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y_dbl, is_newton,
                                         true, scale_step, xtol, ftol, maxfev);
}

TEST_F(algebra_solver_simple_eq_nopara_test, powell) {
  using stan::math::algebra_solver_powell;
  Eigen::VectorXd theta = algebra_solver_powell(simple_eq_functor_nopara(), x,
                                                y_dummy, dat, dummy_dat_int);
  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
}

TEST_F(algebra_solver_non_linear_eq_test, powell_dbl) {
  bool is_newton = false;
  Eigen::VectorXd theta
      = non_linear_eq_test(non_linear_eq_functor(), y_dbl, is_newton);
  EXPECT_FLOAT_EQ(-y_dbl(0), theta(0));
  EXPECT_FLOAT_EQ(-y_dbl(1), theta(1));
  EXPECT_FLOAT_EQ(y_dbl(2), theta(2));
}

TEST_F(algebra_solver_simple_eq_nopara_test, powell_double) {
  using stan::math::algebra_solver_powell;
  Eigen::VectorXd theta = algebra_solver_powell(simple_eq_functor_nopara(), x,
                                                y_dummy, dat, dummy_dat_int);
  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
}

TEST_F(error_message_test, powell_dbl) {
  bool is_newton = false;
  error_conditions_test(non_linear_eq_functor(), y_3, is_newton);
}

TEST(unsolvable_test, powell_dbl) {
  Eigen::VectorXd y(2);
  y << 1, 1;
  unsolvable_test(y);
}

TEST_F(max_steps_test, powell_dbl) {
  bool is_newton = false;
  max_num_steps_test(y, is_newton);
}

TEST_F(degenerate_eq_test, powell_guess1_dbl) {
  using stan::math::algebra_solver_powell;

  // This first initial guess produces the
  // solution x = {8, 8}

  Eigen::VectorXd theta = algebra_solver_powell(degenerate_eq_functor(),
                                                x_guess_1, y_dbl, dat, dat_int);
  EXPECT_FLOAT_EQ(8, theta(0));
  EXPECT_FLOAT_EQ(8, theta(1));
}

TEST_F(degenerate_eq_test, powell_guess2_dbl) {
  using stan::math::algebra_solver_powell;
  // This next initial guess produces the
  // solution x = {5, 5}

  Eigen::VectorXd theta = algebra_solver_powell(degenerate_eq_functor(),
                                                x_guess_2, y_dbl, dat, dat_int);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));
}

// For the next two unit tests,  see if the initial
// guess determines neighborhood of the
// solution, when solutions have different scales,
// using y_scale.

TEST_F(degenerate_eq_test, powell_guess2_scale_dbl) {
  using stan::math::algebra_solver_powell;

  Eigen::VectorXd theta = algebra_solver_powell(
      degenerate_eq_functor(), x_guess_2, y_scale, dat, dat_int);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));
}

TEST_F(degenerate_eq_test, powell_guess_saddle_point_dbl) {
  using stan::math::algebra_solver_powell;

  Eigen::VectorXd theta = algebra_solver_powell(
      degenerate_eq_functor(), x_guess_3, y_scale, dat, dat_int);
  EXPECT_FLOAT_EQ(100, theta(0));
  EXPECT_FLOAT_EQ(100, theta(1));
}

TEST_F(algebra_solver_simple_eq_test, powell) {
  using stan::math::var;
  bool is_newton = false;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y, is_newton);

    std::vector<stan::math::var> y_vec{y(0), y(1), y(2)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++) {
      EXPECT_EQ(J(k, i), g[i]);
    }
  }
}

TEST_F(algebra_solver_simple_eq_test, powell_tuned) {
  using stan::math::var;
  bool is_newton = false;
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

TEST_F(algebra_solver_simple_eq_test, powell_init_is_para) {
  using stan::math::algebra_solver_powell;
  Eigen::VectorXd theta
      = algebra_solver_powell(simple_eq_functor(), x_var, y_dbl, dat, dat_int);
  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
}

TEST_F(algebra_solver_non_linear_eq_test, powell) {
  using stan::math::var;
  bool is_newton = false;
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

TEST_F(error_message_test, powell) {
  using stan::math::var;
  bool is_newton = false;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_2;
  error_conditions_test(non_linear_eq_functor(), y, is_newton);
}

TEST(unsolvable_test, powell) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
  y << 1, 1;
  unsolvable_test(y);
}

TEST_F(max_steps_test, powell) {
  bool is_newton = false;
  max_num_steps_test(y_var, is_newton);
}

TEST_F(degenerate_eq_test, powell_guess1) {
  using stan::math::algebra_solver_powell;
  using stan::math::var;

  // This first initial guess produces the
  // solution x = {8, 8}
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = algebra_solver_powell(
        degenerate_eq_functor(), x_guess_1, y, dat, dat_int);
    EXPECT_FLOAT_EQ(8, theta(0).val());
    EXPECT_FLOAT_EQ(8, theta(1).val());

    std::vector<stan::math::var> y_vec{y(0), y(1)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J1(k, l), g[l], tolerance);
  }
}

TEST_F(degenerate_eq_test, powell_guess2) {
  using stan::math::algebra_solver_powell;
  using stan::math::var;
  // This next initial guess produces the
  // solution x = {5, 5}
  for (int k = 0; k < 1; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = algebra_solver_powell(
        degenerate_eq_functor(), x_guess_2, y, dat, dat_int);
    EXPECT_FLOAT_EQ(5, theta(0).val());
    EXPECT_FLOAT_EQ(5, theta(0).val());

    std::vector<stan::math::var> y_vec{y(0), y(1)};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J2(k, l), g[l], tolerance);
  }
}

TEST_F(variadic_test, powell) {
  using stan::math::var;
  bool is_newton = true;
  for (int k = 0; k < n_x; k++) {
    var y_1 = y_1_dbl;
    var y_2 = y_2_dbl;
    var y_3 = y_3_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = variadic_eq_test(variadic_eq_functor(), A, y_1, y_2, y_3, i,
                           is_newton, scaling_step_size, relative_tolerance,
                           function_tolerance, max_num_steps);

    std::vector<var> y_vec{y_1, y_2, y_3};
    std::vector<double> g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_NEAR(J(k, i), g[i], 1e-6);
  }
}
