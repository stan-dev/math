#include <gtest/gtest.h>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

/* define algebraic functions which get solved */

struct simple_eq_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - y(0) * y(1);
    z(1) = x(1) - y(2);
    return z;
  }
};

struct simple_eq_functor_nopara {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - dat[0] * dat[1];
    z(1) = x(1) - dat[2];
    return z;
  }
};

struct non_linear_eq_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(3);
    z(0) = x(2) - y(2);
    z(1) = x(0) * x(1) - y(0) * y(1);
    z(2) = (x(2) / x(0) + y(2) / y(0));
    return z;
  }
};

struct non_square_eq_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) - y(0);
    z(1) = x(1) * x(2) - y(1);
    return z;
  }
};

struct unsolvable_eq_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
    z(0) = x(0) * x(0) + y(0);
    z(1) = x(1) * x(1) + y(1);
    return z;
  }
};

// Degenerate roots: each solution can either be y(0) or y(1).
struct degenerate_eq_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
    z(0) = (x(0) - y(0)) * (x(1) - y(1));
    z(1) = (x(0) - y(1)) * (x(1) - y(0));
    return z;
  }
};

/* template code for running tests in the prim and rev regime */

template <typename F, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> simple_eq_test(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
    bool tuning = false, double rel_tol = 1e-10, double fun_tol = 1e-6,
    int32_t max_steps = 1e+3) {
  using stan::math::algebra_solver;
  using stan::math::var;

  int n_x = 2;
  Eigen::VectorXd x(n_x);
  x << 1, 1;  // initial guess
  std::vector<double> dummy_dat;
  std::vector<int> dummy_dat_int;

  Eigen::Matrix<T, Eigen::Dynamic, 1> theta;

  if (tuning == false)
    theta = algebra_solver(f, x, y, dummy_dat, dummy_dat_int);
  else
    theta = algebra_solver(f, x, y, dummy_dat, dummy_dat_int, 0, rel_tol,
                           fun_tol, max_steps);

  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));

  return theta;
}

template <typename F, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> non_linear_eq_test(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y) {
  using stan::math::algebra_solver;
  using stan::math::var;

  int n_x = 3;
  Eigen::VectorXd x(n_x);
  x << -3, -3, -3;  // note: need good guess for this one
  std::vector<double> dummy_dat;
  std::vector<int> dummy_dat_int;

  Eigen::Matrix<T, Eigen::Dynamic, 1> theta;

  theta = algebra_solver(f, x, y, dummy_dat, dummy_dat_int);

  return theta;
}

template <typename F, typename T>
inline void error_conditions_test(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y) {
  using stan::math::algebra_solver;

  int n_x = 3;
  Eigen::VectorXd x(n_x);
  x << 1, 1, 1;
  std::vector<double> dat;
  std::vector<int> dat_int;

  std::stringstream err_msg;
  err_msg << "algebra_solver: size of the algebraic system's output (2) "
          << "and size of the vector of unknowns, x, (3) must match in size";
  std::string msg = err_msg.str();
  EXPECT_THROW_MSG(algebra_solver(non_square_eq_functor(), x, y, dat, dat_int),
                   std::invalid_argument, msg);

  Eigen::VectorXd x_bad(static_cast<Eigen::VectorXd::Index>(0));
  std::stringstream err_msg2;
  err_msg2 << "algebra_solver: initial guess has size 0, but "
           << "must have a non-zero size";
  std::string msg2 = err_msg2.str();
  EXPECT_THROW_MSG(algebra_solver(f, x_bad, y, dat, dat_int),
                   std::invalid_argument, msg2);

  double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd x_bad_inf(n_x);
  x_bad_inf << inf, 1, 1;
  EXPECT_THROW_MSG(algebra_solver(f, x_bad_inf, y, dat, dat_int),
                   std::domain_error,
                   "algebra_solver: initial guess is inf, but must "
                   "be finite!");

  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> matrix;
  matrix y_bad_inf(3);
  y_bad_inf << inf, 1, 1;
  EXPECT_THROW_MSG(algebra_solver(f, x, y_bad_inf, dat, dat_int),
                   std::domain_error,
                   "algebra_solver: parameter vector is inf, but must "
                   "be finite!");

  std::vector<double> dat_bad_inf(1);
  dat_bad_inf[0] = inf;
  EXPECT_THROW_MSG(algebra_solver(f, x, y, dat_bad_inf, dat_int),
                   std::domain_error,
                   "algebra_solver: continuous data is inf, but must "
                   "be finite!");

  EXPECT_THROW_MSG(algebra_solver(f, x, y, dat, dat_int, 0, -1, 1e-6, 1e+3),
                   std::invalid_argument, "relative_tolerance");

  EXPECT_THROW_MSG(algebra_solver(f, x, y, dat, dat_int, 0, 1e-6, -1, 1e+3),
                   std::invalid_argument, "function_tolerance");

  EXPECT_THROW_MSG(algebra_solver(f, x, y, dat, dat_int, 0, 1e-6, 1e-6, -1),
                   std::invalid_argument, "max_num_steps");
}

template <typename T>
void inline unsolvable_test(Eigen::Matrix<T, Eigen::Dynamic, 1>& y) {
  using stan::math::algebra_solver;

  Eigen::VectorXd x(2);
  x << 1, 1;
  std::vector<double> dat;
  std::vector<int> dat_int;
  double relative_tolerance = 1e-6, function_tolerance = 1e-6;
  int max_num_steps = 1e+3;

  std::stringstream err_msg;
  err_msg << "algebra_solver: the norm of the algebraic function is: "
          << 1.41421  // sqrt(2)
          << " but should be lower than the function tolerance: "
          << function_tolerance
          << ". Consider decreasing the relative tolerance and increasing"
          << " the max_num_steps.";
  std::string msg = err_msg.str();
  EXPECT_THROW_MSG(
      algebra_solver(unsolvable_eq_functor(), x, y, dat, dat_int, 0,
                     relative_tolerance, function_tolerance, max_num_steps),
      boost::math::evaluation_error, msg);
}

template <typename T>
inline void max_num_steps_test(Eigen::Matrix<T, Eigen::Dynamic, 1>& y) {
  using stan::math::algebra_solver;

  Eigen::VectorXd x(2);
  x << 1, 1;
  std::vector<double> dat;
  std::vector<int> dat_int;

  double relative_tolerance = 1e-6, function_tolerance = 1e-6;
  int max_num_steps = 2;  // very low for test

  std::stringstream err_msg;
  err_msg << "algebra_solver: max number of iterations: " << max_num_steps
          << " exceeded.";
  std::string msg = err_msg.str();
  EXPECT_THROW_MSG(
      algebra_solver(unsolvable_eq_functor(), x, y, dat, dat_int, 0,
                     relative_tolerance, function_tolerance, max_num_steps),
      boost::math::evaluation_error, msg);
}

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> degenerate_test(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& y, const Eigen::VectorXd& x) {
  using stan::math::algebra_solver;

  std::vector<double> dat;
  std::vector<int> dat_int;

  return algebra_solver(degenerate_eq_functor(), x, y, dat, dat_int);
}
