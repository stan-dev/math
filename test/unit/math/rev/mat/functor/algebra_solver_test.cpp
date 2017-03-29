#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <test/unit/util.hpp>


template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
simpleEq(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
          const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
          const std::vector<double>& dat,
          const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
  z(0) = x(0) - y(0) * y(1);
  z(1) = x(1) - y(2);
  return z;
}

template <typename T1, typename T2>
struct simpleEq_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> x_;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> y_;
  std::vector<double> dat_;
  std::vector<int> dat_int_;
  bool x_is_dv_;

public:
  simpleEq_functor() { };

  simpleEq_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                   const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
                   const std::vector<double>& dat,
                   const std::vector<int>& dat_int,
                   const bool& x_is_dv) {
    x_ = x;
    y_ = y;
    dat_ = dat;
    dat_int_ = dat_int;
    x_is_dv_ = x_is_dv;
  }

  template <typename T>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T, T1, T2>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    if (x_is_dv_)
      return simpleEq(x, y_, dat_, dat_int_);
    else
      return simpleEq(x_, x, dat_, dat_int_);
  }
};

TEST(MathMatrix, simple_Eq) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;

  int n_x = 2, n_y = 3;
  for (int k = 0; k < n_x; k++) {
    Eigen::VectorXd x(n_x);
    x << 1, 1;  // initial guess
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
    y << 5, 4, 2;
    std::vector<double> dummy_dat(0);
    std::vector<int> dummy_dat_int(0);

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta;

    theta = algebra_solver(simpleEq_functor<double, double>(),
                           x, y, dummy_dat, dummy_dat_int);

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
}


///////////////////////////////////////////////////////////////////////////////

template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
nonLinearEq(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
            const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
            const std::vector<double>& dat,
            const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(3);
  z(0) = x(2) - y(2);
  z(1) = x(0) * x(1) - y(0) * y(1);
  z(2) = (x(2) / x(0) + y(2) / y(0));
  return z;
}

template <typename T1, typename T2>
struct nonLinearEq_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> x_;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> y_;
  std::vector<double> dat_;
  std::vector<int> dat_int_;
  bool x_is_dv_;

public:
  nonLinearEq_functor() { };

  nonLinearEq_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
                      const std::vector<double>& dat,
                      const std::vector<int>& dat_int,
                      const bool& x_is_dv) {
    x_ = x;
    y_ = y;
    dat_ = dat;
    dat_int_ = dat_int;
    x_is_dv_ = x_is_dv;
  }

  template <typename T>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T, T1, T2>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    if (x_is_dv_)
      return nonLinearEq(x, y_, dat_, dat_int_);
    else
      return nonLinearEq(x_, x, dat_, dat_int_);
  }
};

TEST(MathMatrix, nonLinearEq) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;

  int n_x = 3, n_y = 3;
  for (int k = 0; k < n_x; k++) {
    Eigen::VectorXd x(n_x);
    x << -4, -6, 3;  // initial guess FIX ME: what should these be?
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
    y << 4, 6, 3;
    std::vector<double> dat(1);
    dat[0] = 7;
    std::vector<int> dummy_dat_int(0);

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta;

    theta = algebra_solver(nonLinearEq_functor<double, double>(),
                           x, y, dat, dummy_dat_int);

    EXPECT_FLOAT_EQ(- y(0).val(), theta(0).val());
    EXPECT_FLOAT_EQ(- y(1).val(), theta(1).val());
    EXPECT_FLOAT_EQ(y(2).val(), theta(2).val());

    Eigen::MatrixXd J(n_x, n_y);
    J << -1, 0, 0, 0, -1, 0, 0, 0, 1;

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

///////////////////////////////////////////////////////////////////////////////

template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
nonSquareEq(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
            const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
            const std::vector<double>& dat,
            const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
  z(0) = x(0) - y(0);
  z(1) = x(1) * x(2) - y(1);
  return z;
}

template <typename T1, typename T2>
struct nonSquareEq_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> x_;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> y_;
  std::vector<double> dat_;
  std::vector<int> dat_int_;
  bool x_is_dv_;

public:
  nonSquareEq_functor() { };

  nonSquareEq_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
                      const std::vector<double>& dat,
                      const std::vector<int>& dat_int,
                      const bool& x_is_dv) {
    x_ = x;
    y_ = y;
    dat_ = dat;
    dat_int_ = dat_int;
    x_is_dv_ = x_is_dv;
  }

  template <typename T>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T, T1, T2>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    if (x_is_dv_)
      return nonSquareEq(x, y_, dat_, dat_int_);
    else
      return nonSquareEq(x_, x, dat_, dat_int_);
  }
};

TEST(MathMatrix, error_conditions) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;

  int n_x = 3, n_y = 2;
  Eigen::VectorXd x(n_x);
  x << 1, 1, 1;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
  y << 4, 6;
  std::vector<double> dat(0);
  std::vector<int> dat_int(0);

  std::stringstream err_msg;
  err_msg << "algebra_solver: the ouput of the algebraic system has "
          << "dimension = 2, but should have the same dimension as x "
          << "(the vector of unknowns), which is: 3";
  std::string msg = err_msg.str();
  EXPECT_THROW_MSG(algebra_solver(nonSquareEq_functor<double, double>(),
                                  x, y, dat, dat_int),
                   std::invalid_argument,
                   msg);

  Eigen::VectorXd x_bad(static_cast<Eigen::VectorXd::Index>(0));
  std::stringstream err_msg2;
  err_msg2 << "algebra_solver: initial guess has size 0, but "
            << "must have a non-zero size";
  std::string msg2 = err_msg2.str();
  EXPECT_THROW_MSG(algebra_solver(nonLinearEq_functor<double, double>(),
                                  x_bad, y, dat, dat_int),
                   std::invalid_argument,
                   msg2);

  typedef Eigen::Matrix<var, Eigen::Dynamic, 1> matrix_v;
  matrix_v y_bad(static_cast<matrix_v::Index>(0));
  std::stringstream err_msg3;
  err_msg3 << "algebra_solver: parameter vector has size 0, but "
           << "must have a non-zero size";
  std::string msg3 = err_msg3.str();
  EXPECT_THROW_MSG(algebra_solver(nonLinearEq_functor<double, double>(),
                                  x, y_bad, dat, dat_int),
                   std::invalid_argument,
                   msg3);

  double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd x_bad_inf(n_x);
  x_bad_inf << inf, 1, 1;
  EXPECT_THROW_MSG(algebra_solver(nonLinearEq_functor<double, double>(),
                                  x_bad_inf, y, dat, dat_int),
                   std::domain_error,
                   "algebra_solver: initial guess is inf, but must be finite!");

  matrix_v y_bad_inf(3);
  y_bad_inf << inf, 1, 1;
  EXPECT_THROW_MSG(algebra_solver(nonLinearEq_functor<double, double>(),
                                  x, y_bad_inf, dat, dat_int),
                   std::domain_error,
                   "algebra_solver: parameter vector is inf, but must be finite!");

  std::vector<double> dat_bad_inf(1);
  dat_bad_inf[0] = inf;
  EXPECT_THROW_MSG(algebra_solver(nonLinearEq_functor<double, double>(),
                                  x, y, dat_bad_inf, dat_int),
                   std::domain_error,
                   "algebra_solver: continuous data is inf, but must be finite!");

}
