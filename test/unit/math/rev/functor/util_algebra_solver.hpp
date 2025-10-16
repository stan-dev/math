#include <gtest/gtest.h>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/functor/solve_powell.hpp>
#include <stan/math/rev/functor/solve_newton.hpp>
#include <test/unit/util.hpp>
#include <sstream>
#include <vector>
#include <limits>
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

  void TearDown() override { stan::math::recover_memory(); }

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

  void TearDown() override { stan::math::recover_memory(); }

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

  void TearDown() override { stan::math::recover_memory(); }

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

  void TearDown() override { stan::math::recover_memory(); }

  Eigen::VectorXd y_2;
  Eigen::VectorXd y_3;
};

class max_steps_test : public ::testing::Test {
 protected:
  void SetUp() override {
    y = stan::math::to_vector({1, 1, 1});
    y_var = stan::math::to_vector({1, 1, 1});
  }

  void TearDown() override { stan::math::recover_memory(); }

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

  void TearDown() override { stan::math::recover_memory(); }

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

class variadic_test : public ::testing::Test {
 protected:
  void SetUp() override {
    n_x = 2;
    n_y = 3;

    y_1_dbl = 5;
    y_2_dbl = 4;
    y_3_dbl = 2;

    Eigen::MatrixXd A_(n_x, n_x);
    A_ << 1, 2, 2, 1;
    A = A_;

    x_var = stan::math::to_vector({1, 3});

    Eigen::MatrixXd J_(n_x, n_y);
    J_ << 4, 5, 0, 0, 0, 1;
    J = J_;

    i = 3;
  }

  void TearDown() override { stan::math::recover_memory(); }

  int n_x;
  int n_y;
  Eigen::MatrixXd A_;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> A;
  double y_1_dbl;
  double y_2_dbl;
  double y_3_dbl;
  int i;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x_var;

  double scaling_step_size = 1e-3;
  double relative_tolerance = 1e-6;
  double function_tolerance = 1e-3;
  int max_num_steps = 1e+2;

  Eigen::MatrixXd J;
};

/* wrapper function that either calls the Newton or the Powell solver. */
template <typename F, typename T1, typename T2>
Eigen::Matrix<T2, Eigen::Dynamic, 1> general_algebra_solver(
    bool is_newton, bool use_tol, const F& f,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs = nullptr, double relative_tolerance = 1e-10,
    double scaling_step_size = 1e-3, double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)
  using stan::math::solve_newton;
  using stan::math::solve_newton_tol;
  using stan::math::solve_powell;
  using stan::math::solve_powell_tol;

  Eigen::Matrix<T2, Eigen::Dynamic, 1> theta;

  if (!is_newton) {
    theta = use_tol
                ? solve_powell_tol(f, x, relative_tolerance, function_tolerance,
                                   max_num_steps, msgs, y, dat, dat_int)
                : solve_powell(f, x, msgs, y, dat, dat_int);
  } else {
    theta = use_tol
                ? solve_newton_tol(f, x, scaling_step_size, function_tolerance,
                                   max_num_steps, msgs, y, dat, dat_int)
                : solve_newton(f, x, msgs, y, dat, dat_int);
  }
  return theta;
}

/* wrapper function that either calls the Newton or the Powell solver,
 * using deprecated signatures.
 */
template <typename F, typename T1, typename T2>
Eigen::Matrix<T2, Eigen::Dynamic, 1> general_algebra_solver_non_varia(
    bool is_newton, bool use_tol, const F& f,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs = nullptr, double relative_tolerance = 1e-10,
    double scaling_step_size = 1e-3, double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)
  using stan::math::algebra_solver;
  using stan::math::algebra_solver_newton;

  Eigen::Matrix<T2, Eigen::Dynamic, 1> theta;
  if (!is_newton) {
    theta = algebra_solver(f, x, y, dat, dat_int, msgs, relative_tolerance,
                           function_tolerance, max_num_steps);
  } else {
    theta
        = algebra_solver_newton(f, x, y, dat, dat_int, msgs, scaling_step_size,
                                function_tolerance, max_num_steps);
  }
  return theta;
}

/* define algebraic functions which get solved. */
struct simple_eq_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, std::ostream* pstream__, const T2& y, const T3& dat,
             const T4& dat_int) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - y(0) * y(1);
    z(1) = x(1) - y(2);
    return z;
  }
};

// Same as above, which pstream__ passed as the last argument when using
// an algebra_solver with no variadic arguments.
struct simple_eq_non_varia_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, const T2& y, const T3& dat, const T4& dat_int,
             std::ostream* pstream__) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - y(0) * y(1);
    z(1) = x(1) - y(2);
    return z;
  }
};

struct simple_eq_functor_nopara {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, std::ostream* pstream__, const T2& y, const T3& dat,
             const T4& dat_int) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - dat[0] * dat[1];
    z(1) = x(1) - dat[2];
    return z;
  }
};

struct non_linear_eq_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, std::ostream* pstream__, const T2& y, const T3& dat,
             const T4& dat_int) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(3);
    z(0) = x(2) - y(2);
    z(1) = x(0) * x(1) - y(0) * y(1);
    z(2) = x(2) / x(0) + y(2) / y(0);
    return z;
  }
};

struct non_square_eq_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, std::ostream* pstream__, const T2& y, const T3& dat,
             const T4& dat_int) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - y(0);
    z(1) = x(1) * x(2) - y(1);
    return z;
  }
};

struct unsolvable_eq_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, std::ostream* pstream__, const T2& y, const T3& dat,
             const T4& dat_int) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) * x(0) + y(0);
    z(1) = x(1) * x(1) + y(1);
    return z;
  }
};

// Degenerate roots: each solution can either be y(0) or y(1).
struct degenerate_eq_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T1& x, std::ostream* pstream__, const T2& y, const T3& dat,
             const T4& dat_int) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> z(2);
    z(0) = (x(0) - y(0)) * (x(1) - y(1));
    z(1) = (x(0) - y(1)) * (x(1) - y(0));
    return z;
  }
};

struct variadic_eq_functor {
  template <typename T1, typename T2, typename T3, typename T4>
  inline Eigen::Matrix<stan::return_type_t<T1, T2, T3, T4>, Eigen::Dynamic, 1>
  operator()(const T1& x, std::ostream* pstream__, const T2& A, const T3& y_1,
             const T3& y_2, const T3& y_3, const T4& i) const {
    Eigen::Matrix<stan::return_type_t<T1, T2, T3, T4>, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - y_1 * y_2;
    z(1) = x(1) - y_3;
    for (int j = 0; j < i; ++j)
      z = A * z;
    return z;
  }
};

/* template code for running tests in the prim and rev regime */

template <typename F, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> simple_eq_test(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
    bool is_newton = false, bool tuning = false, double scale_step = 1e-3,
    double rel_tol = 1e-10, double fun_tol = 1e-6, int32_t max_steps = 1e+3) {
  int n_x = 2;
  Eigen::VectorXd x(n_x);
  x << 1, 1;  // initial guess
  std::vector<double> dummy_dat;
  std::vector<int> dummy_dat_int;

  Eigen::Matrix<T, Eigen::Dynamic, 1> theta = general_algebra_solver(
      is_newton, tuning, f, x, y, dummy_dat, dummy_dat_int, &std::cout,
      scale_step, rel_tol, fun_tol, max_steps);

  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));

  return theta;
}

template <typename F, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> simple_eq_non_varia_test(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
    bool is_newton = false, bool tuning = false, double scale_step = 1e-3,
    double rel_tol = 1e-10, double fun_tol = 1e-6, int32_t max_steps = 1e+3) {
  int n_x = 2;
  Eigen::VectorXd x(n_x);
  x << 1, 1;  // initial guess
  std::vector<double> dummy_dat;
  std::vector<int> dummy_dat_int;

  Eigen::Matrix<T, Eigen::Dynamic, 1> theta = general_algebra_solver_non_varia(
      is_newton, tuning, f, x, y, dummy_dat, dummy_dat_int, &std::cout,
      scale_step, rel_tol, fun_tol, max_steps);

  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));

  return theta;
}

template <typename F, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> non_linear_eq_test(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
    int solver_type = 0, bool use_tol = false) {
  int n_x = 3;
  Eigen::VectorXd x(n_x);
  x << -3, -3, -3;  // note: need good guess for this one
  std::vector<double> dummy_dat;
  std::vector<int> dummy_dat_int;

  Eigen::Matrix<T, Eigen::Dynamic, 1> theta;

  theta = general_algebra_solver(solver_type, use_tol, f, x, y, dummy_dat,
                                 dummy_dat_int);

  return theta;
}

template <typename F, typename T>
inline void error_conditions_test(const F& f,
                                  const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                                  int solver_type = 0, bool use_tol = true) {
  int n_x = 3;
  Eigen::VectorXd x(n_x);
  x << 1, 1, 1;
  std::vector<double> dat;
  std::vector<int> dat_int;

  std::stringstream err_msg;
  err_msg << "size of the algebraic system's output (2) "
          << "and size of the vector of unknowns, x, (3) must match in size";
  std::string msg = err_msg.str();
  EXPECT_THROW_MSG(
      general_algebra_solver(solver_type, use_tol, non_square_eq_functor(), x,
                             y, dat, dat_int),
      std::invalid_argument, msg);

  Eigen::VectorXd x_bad(static_cast<Eigen::VectorXd::Index>(0));
  std::stringstream err_msg2;
  err_msg2 << "initial guess has size 0, but "
           << "must have a non-zero size";
  std::string msg2 = err_msg2.str();
  EXPECT_THROW_MSG(
      general_algebra_solver(solver_type, use_tol, f, x_bad, y, dat, dat_int),
      std::invalid_argument, msg2);

  double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd x_bad_inf(n_x);
  x_bad_inf << inf, 1, 1;
  EXPECT_THROW_MSG(general_algebra_solver(solver_type, use_tol, f, x_bad_inf, y,
                                          dat, dat_int),
                   std::domain_error,
                   "initial guess[1] is inf, but must be finite!");

  if (solver_type == 0) {
    EXPECT_THROW_MSG(general_algebra_solver(solver_type, use_tol, f, x, y, dat,
                                            dat_int, 0, -1, 1e-3, 1e-6, 1e+3),
                     std::domain_error, "relative_tolerance");
  }
  if (solver_type == 1) {
    EXPECT_THROW_MSG(general_algebra_solver(solver_type, use_tol, f, x, y, dat,
                                            dat_int, 0, 1e-10, -1, 1e-6, 1e+3),
                     std::domain_error, "scaling_step_size");
  }
  EXPECT_THROW_MSG(general_algebra_solver(solver_type, use_tol, f, x, y, dat,
                                          dat_int, 0, 1e-10, 1e-3, -1, 1e+3),
                   std::domain_error, "function_tolerance");
  EXPECT_THROW_MSG(general_algebra_solver(solver_type, use_tol, f, x, y, dat,
                                          dat_int, 0, 1e-10, 1e-3, 1e-6, -1),
                   std::domain_error, "max_num_steps");
}

template <typename T>
void inline unsolvable_test(Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                            int solver_type = 0, bool use_tol = false) {
  Eigen::VectorXd x(2);
  x << 1, 1;
  std::vector<double> dat;
  std::vector<int> dat_int;
  double relative_tolerance = 1e-6, function_tolerance = 1e-6,
         scaling_step_size = 1e-3;
  int max_num_steps = 1e+3;

  std::stringstream err_msg;
  err_msg << "algebra_solver: the norm of the algebraic function is "
          << 1.41421  // sqrt(2)
          << " but should be lower than the function tolerance: "
          << function_tolerance
          << ". Consider decreasing the relative tolerance and increasing"
          << " max_num_steps.";
  std::string msg = err_msg.str();
  EXPECT_THROW_MSG(general_algebra_solver(
                       solver_type, use_tol, unsolvable_eq_functor(), x, y, dat,
                       dat_int, 0, relative_tolerance, scaling_step_size,
                       function_tolerance, max_num_steps),
                   std::domain_error, msg);
}

template <typename T>
void inline unsolvable_flag_test(Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                                 int solver_type = 0, bool use_tol = false) {
  Eigen::VectorXd x(2);
  x << 1, 1;
  std::vector<double> dat;
  std::vector<int> dat_int;
  std::stringstream err_msg;
  err_msg << "The linear solver’s setup function failed in an unrecoverable "
             "manner.";  // NOLINT
  std::string msg = err_msg.str();
  EXPECT_THROW_MSG(
      general_algebra_solver(solver_type, use_tol, unsolvable_eq_functor(), x,
                             y, dat, dat_int),
      std::runtime_error, msg);
}

template <typename T>
inline void max_num_steps_test(Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                               int solver_type = 0, bool use_tol = false) {
  Eigen::VectorXd x(3);
  x << 1, 1, 1;
  std::vector<double> dat;
  std::vector<int> dat_int;

  double relative_tolerance = 1e-6, function_tolerance = 1e-6,
         scaling_step = 1e-3;
  int max_num_steps = 2;  // very low for test

  EXPECT_THROW(
      general_algebra_solver(
          solver_type, use_tol, non_linear_eq_functor(), x, y, dat, dat_int, 0,
          scaling_step, relative_tolerance, function_tolerance, max_num_steps),
      std::domain_error);
}

Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> variadic_eq_impl_test(
    const Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> A,
    const stan::math::var& y_1, const stan::math::var& y_2,
    const stan::math::var& y_3, const int i, bool is_newton = false,
    bool use_tol = false, double scaling_step_size = 1e-3,
    double relative_tolerance = 1e-10, double function_tolerance = 1e-6,
    int32_t max_num_steps = 1e+3) {
  using stan::math::solve_newton;
  using stan::math::solve_powell;
  using stan::math::solve_powell_tol;
  using stan::math::var;

  Eigen::VectorXd x(2);
  x << 1, 1;  // initial guess

  Eigen::Matrix<var, Eigen::Dynamic, 1> theta;

  theta = is_newton
              ? use_tol ? solve_newton_tol(variadic_eq_functor(), x,
                                           relative_tolerance,
                                           function_tolerance, max_num_steps,
                                           &std::cout, A, y_1, y_2, y_3, i)
                        : solve_newton(variadic_eq_functor(), x, &std::cout, A,
                                       y_1, y_2, y_3, i)
              : use_tol ? solve_powell_tol(variadic_eq_functor(), x,
                                           relative_tolerance,
                                           function_tolerance, max_num_steps,
                                           &std::cout, A, y_1, y_2, y_3, i)
                        : solve_powell(variadic_eq_functor(), x, &std::cout, A,
                                       y_1, y_2, y_3, i);

  EXPECT_NEAR(20, value_of(theta(0)), 1e-6);
  EXPECT_NEAR(2, value_of(theta(1)), 1e-6);

  return theta;
}
