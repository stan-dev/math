#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_solver_powell.hpp>
#include <stan/math/rev/functor/algebra_solver_newton.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/functor/util_algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Every test exists in four versions for the cases
// where y (the auxiliary parameters) are passed as
// data (double type) or parameters (var types),
// and the cases where the solver is based on Powell's
// or Newton's method.

class algebra_solver_simple_eq_test : public ::testing::Test {
 protected:
  void SetUp() override {
    n_x = 2;
    n_y = 3;

    y_dbl = stan::math::to_vector({5, 4, 2});
    Eigen::MatrixXd J_(n_x, n_y);
    J_ << 4, 5, 0, 0, 0, 1;
    J = J_;

    x_var = stan::math::to_vector({1, 1});
  }

  int n_x;
  int n_y;
  Eigen::VectorXd y_dbl;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x_var;
  std::vector<double> dat;
  std::vector<int> dat_int;
  double scale_step = 1e-3;
  double xtol = 1e-6;
  double ftol = 1e-3;
  int maxfev = 1e+2;

  Eigen::MatrixXd J;
};

class algebra_solver_simple_eq_nopara_test : public ::testing::Test {
 protected:
  void SetUp() override { x = stan::math::to_vector({1, 1}); }

  int n_x = 2;
  Eigen::VectorXd x;
  std::vector<double> dat = {5, 4, 2};
  Eigen::VectorXd y_dummy;
  std::vector<int> dummy_dat_int;
};

class algebra_solver_non_linear_eq_test : public ::testing::Test {
 protected:
  void SetUp() override {
    y_dbl = stan::math::to_vector({4, 6, 3});
    Eigen::MatrixXd J_(n_x, n_y);
    J_ << -1, 0, 0, 0, -1, 0, 0, 0, 1;
    J = J_;
  }
  int n_x = 3;
  int n_y = 3;
  double err = 1e-11;

  Eigen::VectorXd y_dbl;
  Eigen::MatrixXd J;
};

class error_message_test : public ::testing::Test {
 protected:
  void SetUp() override {
    y_2 = stan::math::to_vector({4, 6});
    y_3 = stan::math::to_vector({4, 6, 3});
  }

  Eigen::VectorXd y_2;
  Eigen::VectorXd y_3;
};

class max_steps_test : public ::testing::Test {
 protected:
  void SetUp() override {
    y = stan::math::to_vector({1, 1, 1});
    y_var = stan::math::to_vector({1, 1, 1});
  }

  Eigen::VectorXd y;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y_var;
};

class degenerate_eq_test : public ::testing::Test {
 protected:
  void SetUp() override {
    using stan::math::to_vector;
    y_dbl = to_vector({5, 8});
    y_scale = to_vector({5, 100});
    x_guess_1 = to_vector({10, 1});
    x_guess_2 = to_vector({1, 1});
    x_guess_3 = to_vector({5, 100});  // 80, 80

    Eigen::MatrixXd J_(n_x, n_y);
    J_ << 0, 1, 0, 1;
    J1 = J_;
    J_ << 1, 0, 1, 0;
    J2 = J_;
  }

  int n_x = 2;
  int n_y = 2;
  double tolerance = 1e-10;
  Eigen::VectorXd y_dbl;
  Eigen::VectorXd y_scale;
  Eigen::VectorXd x_guess_1;
  Eigen::VectorXd x_guess_2;
  Eigen::VectorXd x_guess_3;
  Eigen::MatrixXd J1;
  Eigen::MatrixXd J2;
  std::vector<double> dat;
  std::vector<int> dat_int;
};

//////////////////////////////////////////////////////////////////////////
// Tests for powell solver.

TEST_F(algebra_solver_simple_eq_test, powell) {
  using stan::math::var;
  bool is_newton = false;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y, is_newton);

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
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

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

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

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_NEAR(J(k, i), g[i], err);
  }
}

TEST_F(algebra_solver_non_linear_eq_test, powell_dbl) {
  bool is_newton = false;
  Eigen::VectorXd theta
      = non_linear_eq_test(non_linear_eq_functor(), y_dbl, is_newton);
  EXPECT_FLOAT_EQ(-y_dbl(0), theta(0));
  EXPECT_FLOAT_EQ(-y_dbl(1), theta(1));
  EXPECT_FLOAT_EQ(y_dbl(2), theta(2));
}

TEST_F(error_message_test, powell) {
  using stan::math::var;
  bool is_newton = false;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_2;
  error_conditions_test(non_linear_eq_functor(), y, is_newton);
}

TEST_F(error_message_test, powell_dbl) {
  bool is_newton = false;
  error_conditions_test(non_linear_eq_functor(), y_3, is_newton);
}

TEST(unsolvable_test, powell) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(2);
  y << 1, 1;
  unsolvable_test(y);
}

TEST(unsolvable_test, powell_dbl) {
  Eigen::VectorXd y(2);
  y << 1, 1;
  unsolvable_test(y);
}

TEST_F(max_steps_test, powell) {
  bool is_newton = false;
  max_num_steps_test(y_var, is_newton);
}

TEST_F(max_steps_test, powell_dbl) {
  bool is_newton = false;
  max_num_steps_test(y, is_newton);
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

    AVEC y_vec = createAVEC(y(0), y(1));
    VEC g;
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

    AVEC y_vec = createAVEC(y(0), y(1));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J2(k, l), g[l], tolerance);
  }
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

// unit test to demo issue #696
// system functor init bug issue #696
TEST(MathMatrixRevMat, system_functor_constructor) {
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
//////////////////////////////////////////////////////////////////////////
// Tests for newton solver.

TEST_F(algebra_solver_simple_eq_test, newton) {
  using stan::math::var;
  bool is_newton = true;
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
        = simple_eq_test(simple_eq_functor(), y, is_newton);

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
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

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

TEST_F(algebra_solver_simple_eq_test, newton_dbl) {
  bool is_newton = true;
  Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y_dbl, is_newton);
}

TEST_F(algebra_solver_simple_eq_test, newton_tuned_dbl) {
  bool is_newton = true;
  Eigen::VectorXd theta = simple_eq_test(simple_eq_functor(), y_dbl, is_newton,
                                         true, scale_step, xtol, ftol, maxfev);
}

TEST_F(algebra_solver_simple_eq_nopara_test, newton) {
  using stan::math::algebra_solver_newton;
  Eigen::VectorXd theta = algebra_solver_newton(simple_eq_functor_nopara(), x,
                                                y_dummy, dat, dummy_dat_int);
  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));
}

TEST_F(algebra_solver_simple_eq_test, newton_init_is_para) {
  using stan::math::algebra_solver_newton;
  Eigen::VectorXd theta
      = algebra_solver_newton(simple_eq_functor(), x_var, y_dbl, dat, dat_int);
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

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_NEAR(J(k, i), g[i], err);
  }
}

TEST_F(algebra_solver_non_linear_eq_test, newton_dbl) {
  bool is_newton = true;
  Eigen::VectorXd theta
      = non_linear_eq_test(non_linear_eq_functor(), y_dbl, is_newton);
  EXPECT_FLOAT_EQ(-y_dbl(0), theta(0));
  EXPECT_FLOAT_EQ(-y_dbl(1), theta(1));
  EXPECT_FLOAT_EQ(y_dbl(2), theta(2));
}

TEST_F(error_message_test, newton) {
  using stan::math::var;
  bool is_newton = true;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_2;
  error_conditions_test(non_linear_eq_functor(), y, is_newton);
}

TEST_F(error_message_test, newton_dbl) {
  bool is_newton = true;
  error_conditions_test(non_linear_eq_functor(), y_3, is_newton);
}

TEST_F(max_steps_test, newton) {
  bool is_newton = true;
  max_num_steps_test(y_var, is_newton);
}

TEST_F(max_steps_test, newton_dbl) {
  bool is_newton = true;
  max_num_steps_test(y, is_newton);
}

TEST(MathMatrixRevMat, unsolvable_flag_newton) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y(2);
  y << 1, 1;

  unsolvable_flag_test(y);
}

TEST(MathMatrixRevMat, unsolvable_flag_newton_dbl) {
  Eigen::VectorXd y(2);
  y << 1, 1;

  unsolvable_flag_test(y);
}

TEST_F(degenerate_eq_test, newton_guess1) {
  using stan::math::algebra_solver_newton;
  // using stan::math::sum;
  using stan::math::var;

  // This first initial guess produces the
  // solution x = {8, 8}
  for (int k = 0; k < n_x; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = algebra_solver_newton(
        degenerate_eq_functor(), x_guess_1, y, dat, dat_int);
    EXPECT_FLOAT_EQ(8, theta(0).val());
    EXPECT_FLOAT_EQ(8, theta(1).val());

    AVEC y_vec = createAVEC(y(0), y(1));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J1(k, l), g[l], tolerance);
  }
}

TEST_F(degenerate_eq_test, newton_guess2) {
  using stan::math::algebra_solver_newton;
  using stan::math::var;
  // This next initial guess produces the
  // solution x = {5, 5}
  for (int k = 0; k < 1; k++) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y = y_dbl;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta = algebra_solver_newton(
        degenerate_eq_functor(), x_guess_2, y, dat, dat_int);
    EXPECT_FLOAT_EQ(5, theta(0).val());
    EXPECT_FLOAT_EQ(5, theta(0).val());

    AVEC y_vec = createAVEC(y(0), y(1));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int l = 0; l < n_y; l++)
      EXPECT_NEAR(J2(k, l), g[l], tolerance);
  }
}

TEST_F(degenerate_eq_test, newton_guess1_dbl) {
  using stan::math::algebra_solver_newton;

  // This first initial guess produces the
  // solution x = {8, 8}

  Eigen::VectorXd theta = algebra_solver_newton(degenerate_eq_functor(),
                                                x_guess_1, y_dbl, dat, dat_int);
  EXPECT_FLOAT_EQ(8, theta(0));
  EXPECT_FLOAT_EQ(8, theta(1));
}

TEST_F(degenerate_eq_test, newton_guess2_dbl) {
  using stan::math::algebra_solver_newton;
  // This next initial guess produces the
  // solution x = {5, 5}

  Eigen::VectorXd theta = algebra_solver_newton(degenerate_eq_functor(),
                                                x_guess_2, y_dbl, dat, dat_int);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));
}

// For the next two unit tests,  see if the initial
// guess determines neighborhood of the
// solution, when solutions have different scales,
// using y_scale.

TEST_F(degenerate_eq_test, newton_guess2_scale_dbl) {
  using stan::math::algebra_solver_newton;

  Eigen::VectorXd theta = algebra_solver_newton(
      degenerate_eq_functor(), x_guess_2, y_scale, dat, dat_int);
  EXPECT_FLOAT_EQ(5, theta(0));
  EXPECT_FLOAT_EQ(5, theta(1));
}

TEST_F(degenerate_eq_test, newton_guess_saddle_point_dbl) {
  // Newton solver fails this test because the initial point is
  // a saddle point.
  using stan::math::algebra_solver_newton;
  std::stringstream err_msg;
  err_msg << "algebra_solver failed with error flag -11.";
  std::string msg = err_msg.str();

  EXPECT_THROW_MSG(algebra_solver_newton(degenerate_eq_functor(), x_guess_3,
                                         y_scale, dat, dat_int),
                   std::runtime_error, msg);
}
