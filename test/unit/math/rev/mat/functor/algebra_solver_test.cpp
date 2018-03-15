#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
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

TEST(MathMatrix, simple_Eq) {
  using stan::math::var;

  int n_x = 2, n_y = 3;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
    y << 5, 4, 2;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y);

    Eigen::MatrixXd J(n_x, n_y);
    J << 4, 5, 0, 0, 0, 1;

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

TEST(MathMatrix, simple_Eq_dbl) {
  Eigen::VectorXd y(3);
  y << 5, 4, 2;
  Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y);
}

TEST(MathMatrix, simple_Eq_tuned) {
  using stan::math::var;

  int n_x = 2, n_y = 3;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
    y << 5, 4, 2;
    double xtol = 1e-6, ftol = 1e-6;
    int maxfev = 1e+4;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y, true, xtol, ftol, maxfev);

    Eigen::MatrixXd J(n_x, n_y);
    J << 4, 5, 0, 0, 0, 1;

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

TEST(MathMatrix, simple_Eq_tuned_dbl) {
  int n_y = 3;
  Eigen::VectorXd y(n_y);
  y << 5, 4, 2;
  double xtol = 1e-6, ftol = 1e-6;
  int maxfev = 1e+4;

  Eigen::VectorXd theta
      = simple_eq_test(simple_eq_functor(), y, true, xtol, ftol, maxfev);
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

  Eigen::Matrix<double, Eigen::Dynamic, 1> theta;

  theta = stan::math::algebra_solver(simple_eq_functor_nopara(), x, y_dummy,
                                     dat, dummy_dat_int);

  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
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

  Eigen::VectorXd theta
      = stan::math::algebra_solver(simple_eq_functor(), x, y, dat, dat_int);

  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
}

TEST(MathMatrix, non_linear_eq) {
  using stan::math::var;

  int n_x = 3, n_y = 3;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
    y << 4, 6, 3;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = non_linear_eq_test(non_linear_eq_functor(), y);

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

TEST(MathMatrix, nonLinearEq_dbl) {
  int n_y = 3;
  Eigen::VectorXd y(n_y);
  y << 4, 6, 3;

  Eigen::VectorXd theta;
  theta = non_linear_eq_test(non_linear_eq_functor(), y);

  EXPECT_FLOAT_EQ(-y(0), theta(0));
  EXPECT_FLOAT_EQ(-y(1), theta(1));
  EXPECT_FLOAT_EQ(y(2), theta(2));
}

TEST(MathMatrix, error_conditions) {
  using stan::math::var;

  int n_y = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
  y << 4, 6;

  error_conditions_test(non_linear_eq_functor(), y);
}

TEST(MathMatrix, error_conditions_dbl) {
  int n_y = 2;
  Eigen::VectorXd y(n_y);
  y << 4, 6;

  error_conditions_test(non_linear_eq_functor(), y);
}

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

TEST(MathMatrix, degenerate) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;

  int n_x = 2, n_y = 2;

  // This first initial guess produces the
  // solution x = {8, 8}
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
    y << 5, 8;
    Eigen::VectorXd x(2);
    x << 10, 1;  // Initial Guess
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = degenerate_test(y, x);
    EXPECT_EQ(8, theta(0));
    EXPECT_EQ(8, theta(1));

    Eigen::MatrixXd J(n_x, n_y);
    J << 0, 1, 0, 1;

    AVEC y_vec = createAVEC(y(0), y(1));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_EQ(J(k, l), g[l]);
  }

  // This next initial guess produces the
  // solution x = {5, 5}
  for (int k = 0; k < 1; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
    y << 5, 8;
    Eigen::VectorXd x(2);
    x << 1, 1;  // Initial Guess
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = degenerate_test(y, x);
    EXPECT_FLOAT_EQ(5, theta(0).val());
    EXPECT_FLOAT_EQ(5, theta(0).val());

    Eigen::MatrixXd J(n_x, n_y);
    J << 1, 0, 1, 0;

    AVEC y_vec = createAVEC(y(0), y(1));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J(k, l), g[l], 1e-13);
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
  theta = degenerate_test(y, x);
  EXPECT_EQ(8, theta(0));
  EXPECT_EQ(8, theta(1));

  // This next initial guess produces the
  // solution x = {5, 5}
  x << 1, 1;  // Initial Guess
  theta = degenerate_test(y, x);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));

  // See if the initial guess determines neighborhood of the
  // solution, when solutions have different scales.
  y << 5, 100;
  x << 1, 1;  // initial guess
  theta = degenerate_test(y, x);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));

  x << 120, 120;  // Initial guess
  theta = degenerate_test(y, x);
  EXPECT_FLOAT_EQ(100, theta(0));
  EXPECT_FLOAT_EQ(100, theta(1));
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
