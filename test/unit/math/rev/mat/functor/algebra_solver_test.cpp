#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_newton.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/functor/util_algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

// Every test exists in duplicate to test the case
// where y (the auxiliary parameters) are passed as
// data (double type) or parameters (var types).
// Within each of these tests, the dogleg and the Newton
// solver are tested.

TEST(MathMatrix, simple_Eq) {
  using stan::math::var;

  int n_x = 2, n_y = 3;
  bool is_newton;

  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    for (int k = 0; k < n_x; k++) {
      Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
      y << 5, 4, 2;

      Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y, is_newton);

      Eigen::MatrixXd J(n_x, n_y);
      J << 4, 5, 0, 0, 0, 1;

      AVEC y_vec = createAVEC(y(0), y(1), y(2));
      VEC g;
      theta(k).grad(y_vec, g);

      for (int i = 0; i < n_y; i++) EXPECT_EQ(J(k, i), g[i]);
    }
  }
}

TEST(MathMatrix, simple_Eq_dbl) {
  Eigen::VectorXd y(3);
  y << 5, 4, 2;
  for (int is_newton = 0; is_newton <= 1; is_newton++)
    Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y, is_newton);
}

TEST(MathMatrix, simple_Eq_tuned) {
  using stan::math::var;

  int n_x = 2, n_y = 3;
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    for (int k = 0; k < n_x; k++) {
      Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
      y << 5, 4, 2;
      double xtol = 1e-6, ftol = 1e-6;
      int maxfev = 1e+4;

      Eigen::Matrix<var, Eigen::Dynamic, 1> theta
          = simple_eq_test(simple_eq_functor(), y, is_newton,
                           true, xtol, ftol, maxfev);

      Eigen::MatrixXd J(n_x, n_y);
      J << 4, 5, 0, 0, 0, 1;

      AVEC y_vec = createAVEC(y(0), y(1), y(2));
      VEC g;
      theta(k).grad(y_vec, g);

      for (int i = 0; i < n_y; i++)
        EXPECT_EQ(J(k, i), g[i]);
    }
  }
}

TEST(MathMatrix, simple_Eq_tuned_dbl) {
  int n_y = 3;
  Eigen::VectorXd y(n_y);
  y << 5, 4, 2;
  double xtol = 1e-6, ftol = 1e-6;
  int maxfev = 1e+4;

  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    Eigen::VectorXd theta
        = simple_eq_test(simple_eq_functor(), y, is_newton,
                         true, xtol, ftol, maxfev);
  }
}

TEST(MathMatrix, simple_Eq_nopara) {
  std::vector<double> dat(3);
  dat[0] = 5;
  dat[1] = 4;
  dat[2] = 2;

  int n_x = 2;
  Eigen::VectorXd x(n_x);
  x << 1, 1;  // initial guess
  Eigen::VectorXd y_dummy;
  std::vector<int> dummy_dat_int;

  // Eigen::Matrix<double, Eigen::Dynamic, 1> theta;
  Eigen::VectorXd theta;

  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    theta = general_algebra_solver(is_newton, simple_eq_functor_nopara(), x,
                                   y_dummy, dat, dummy_dat_int);
    EXPECT_EQ(20, theta(0));
    EXPECT_EQ(2, theta(1));
  }
}

TEST(MathMatrix, simple_Eq_init_is_para) {
  using stan::math::var;
  Eigen::VectorXd y(3);
  y << 5, 4, 2;

  int n_x = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> x(n_x);
  x << 1, 1;

  std::vector<double> dat;
  std::vector<int> dat_int;

  for (int is_newton = 0; is_newton < 1; is_newton++) {
  Eigen::VectorXd theta
      = general_algebra_solver(is_newton, simple_eq_functor(), x, y,
                               dat, dat_int);
    EXPECT_EQ(20, theta(0));
    EXPECT_EQ(2, theta(1));
  }
}

TEST(MathMatrix, non_linear_eq) {
  using stan::math::var;

  int n_x = 3, n_y = 3;
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    for (int k = 0; k < n_x; k++) {
      Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
      y << 4, 6, 3;

      Eigen::Matrix<var, Eigen::Dynamic, 1> theta
          = non_linear_eq_test(non_linear_eq_functor(), y, is_newton);

      EXPECT_FLOAT_EQ(-y(0).val(), theta(0).val());
      EXPECT_FLOAT_EQ(-y(1).val(), theta(1).val());
      EXPECT_FLOAT_EQ(y(2).val(), theta(2).val());

      Eigen::MatrixXd J(n_x, n_y);
      J << -1, 0, 0, 0, -1, 0, 0, 0, 1;

      AVEC y_vec = createAVEC(y(0), y(1), y(2));
      VEC g;
      theta(k).grad(y_vec, g);

      double err = 1e-11;
      for (int i = 0; i < n_y; i++)
        EXPECT_NEAR(J(k, i), g[i], err);
    }
  }
}

TEST(MathMatrix, nonLinearEq_dbl) {
  int n_y = 3;
  Eigen::VectorXd y(n_y);
  y << 4, 6, 3;

  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    Eigen::VectorXd 
      theta = non_linear_eq_test(non_linear_eq_functor(), y, is_newton);

    EXPECT_FLOAT_EQ(-y(0), theta(0));
    EXPECT_FLOAT_EQ(-y(1), theta(1));
    EXPECT_FLOAT_EQ(y(2), theta(2));
  }
}

TEST(MathMatrix, error_conditions) {
  using stan::math::var;

  int n_y = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
  y << 4, 6;
  for (int is_newton = 0; is_newton <= 1; is_newton ++)
    error_conditions_test(non_linear_eq_functor(), y, is_newton);
}

TEST(MathMatrix, error_conditions_dbl) {
  int n_y = 2;
  Eigen::VectorXd y(n_y);
  y << 4, 6;
  for (int is_newton = 0; is_newton <= 1; is_newton++)
    error_conditions_test(non_linear_eq_functor(), y, is_newton);
}
/*
TEST(MathMatrix, unsolvable) {
  using stan::math::var;

  Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
  y << 1, 1;  // should be positive

  unsolvable_test(y);
}

TEST(MathMatrix, unsolvable_dbl) {
  Eigen::VectorXd y(2);
  y << 1, 1;

  unsolvable_test(y);
}

TEST(MathMatrix, max_num_steps) {
  using stan::math::var;

  Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
  y << 1, 1;
  max_num_steps_test(y);
}

TEST(MathMatrix, max_num_steps_dbl) {
  Eigen::VectorXd y(2);
  y << 1, 1;
  max_num_steps_test(y);
}
*/
TEST(MathMatrix, degenerate) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;

  int n_x = 2, n_y = 2;
  double tolerance = 1e-10;

  // This first initial guess produces the
  // solution x = {8, 8}
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    for (int k = 0; k < n_x; k++) {
      Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
      y << 5, 8;
      Eigen::VectorXd x(2);
      x << 10, 1;  // Initial Guess
      Eigen::Matrix<var, Eigen::Dynamic, 1> 
        theta = degenerate_test(y, x, is_newton);
      EXPECT_FLOAT_EQ(8, theta(0).val());
      EXPECT_FLOAT_EQ(8, theta(1).val());

      Eigen::MatrixXd J(n_x, n_y);
      J << 0, 1, 0, 1;

      AVEC y_vec = createAVEC(y(0), y(1));
      VEC g;
      theta(k).grad(y_vec, g);

      for (int l = 0; l < n_y; l++)
        EXPECT_NEAR(J(k, l), g[l], tolerance);
    }
  }

  // This next initial guess produces the
  // solution x = {5, 5}
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    for (int k = 0; k < 1; k++) {
      Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
      y << 5, 8;
      Eigen::VectorXd x(2);
      x << 1, 1;  // Initial Guess
      Eigen::Matrix<var, Eigen::Dynamic, 1> 
        theta = degenerate_test(y, x, is_newton);
      EXPECT_FLOAT_EQ(5, theta(0).val());
      EXPECT_FLOAT_EQ(5, theta(0).val());

    Eigen::MatrixXd J(n_x, n_y);
    J << 1, 0, 1, 0;

      AVEC y_vec = createAVEC(y(0), y(1));
      VEC g;
      theta(k).grad(y_vec, g);

      for (int l = 0; l < n_y; l++)
        EXPECT_NEAR(J(k, l), g[l], tolerance);
    }
  }
}

TEST(MathMatrix, degenerate_dbl) {
  using stan::math::algebra_solver;
  using stan::math::var;

  // This first initial guess produces the
  // solution x = {8, 8}
  Eigen::VectorXd y(2);
  y << 5, 8;
  Eigen::VectorXd x(2);
  x << 10, 1;  // Initial Guess
  Eigen::VectorXd theta;
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    theta = degenerate_test(y, x, is_newton);
    EXPECT_FLOAT_EQ(8, theta(0));
    EXPECT_FLOAT_EQ(8, theta(1));
  }

  // This next initial guess produces the
  // solution x = {5, 5}
  x << 1, 1;  // Initial Guess
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    theta = degenerate_test(y, x, is_newton);
    EXPECT_FLOAT_EQ(5, theta(0));
    EXPECT_FLOAT_EQ(5, theta(1));
  }

  // See if the initial guess determines neighborhood of the
  // solution, when solutions have different scales.
  y << 5, 100;
  x << 1, 1;  // initial guess
  theta = degenerate_test(y, x);
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    EXPECT_FLOAT_EQ(5, theta(0));
    EXPECT_FLOAT_EQ(5, theta(1));
  }

  x << 120, 120;  // Initial guess
  for (int is_newton = 0; is_newton <= 1; is_newton++) {
    theta = degenerate_test(y, x);
    EXPECT_FLOAT_EQ(100, theta(0));
    EXPECT_FLOAT_EQ(100, theta(1));
  }
}

// unit test to demo issue #696
// system functor init bug issue #696
TEST(MathMatrix, system_functor_constructor) {
  using stan::math::system_functor;

  Eigen::VectorXd y(2);
  y << 5, 8;
  Eigen::VectorXd x(2);
  x << 10, 1;
  std::vector<double> dat{0.0, 0.0};
  std::vector<int> dat_int{0, 0};
  std::ostream* msgs = 0;
  int f = 99;

  system_functor<int, double, double, true> fs(f, x, y, dat, dat_int, msgs);

  EXPECT_EQ(fs.f_, f);
}

///////////////////////////////////////////////////////////////////////////////
// Basic tests for Newton solver
/*TEST(MathMatrix, simple_Eq_newton) {
  using stan::math::var;
  using stan::math::algebra_solver_newton;

  int n_x = 2, n_y = 3;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
    y << 5, 4, 2;
    int n_x = 2;
    Eigen::VectorXd x(n_x);
    x << 1, 1;  // initial guess
    std::vector<double> dummy_dat;
    std::vector<int> dummy_dat_int;
    double rel_tol = 1e-6, fun_tol = 1e-6;
    int max_steps = 1e+4;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta;

    // TODO - test without tuning parameters (note: rel_tol is
    // a dummy parameter, and not used in the algorithm)
    theta = algebra_solver_newton(simple_eq_functor(), x, y,
                                  dummy_dat, dummy_dat_int, 0,
                                   rel_tol, fun_tol, max_steps);

    EXPECT_EQ(20, theta(0));
    EXPECT_EQ(2, theta(1));

    Eigen::MatrixXd J(n_x, n_y);
    J << 4, 5, 0, 0, 0, 1;

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}  */
