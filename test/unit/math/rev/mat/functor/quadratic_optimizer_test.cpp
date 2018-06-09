#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/quadratic_optimizer.hpp>
#include <stan/math/rev/mat/functor/eiquadprog.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

using stan::math::quadratic_optimizer;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Dynamic;
using std::cout;

TEST(MathMatrix, eiquadprog) {
// implement example from 
// http://www.labri.fr/perso/guenneba/code/QuadProg/example.cpp.

  using Eigen::solve_quadprog;
  
  MatrixXd G(3, 3);
  G << 2.1, 0.0, 1.0,
       1.5, 2.2, 0.0,
       1.2, 1.3, 3.1;
  
  VectorXd g0(3);
  g0 << 6, 1, 1;
  
  MatrixXd CE(3, 1);
  CE << 1, 2, -1;
  
  VectorXd ce0(1);
  ce0(0) = -4;
  
  MatrixXd CI(3, 4);
  CI << 1, 0, 0, -1,
        0, 1, 0, -1,
        0, 0, 1,  0;

  VectorXd ci0(4);
  ci0 << 0, 0, 0, 10;
  // VectorXd ci0(3);
  // ci0 << -2, -2, -2;  // test inequality constraint
  
  // MatrixXd CI = Eigen::MatrixXd::Identity(3, 3);
  // VectorXd ci0(3);
  // ci0 << -1, -1, -1;

  VectorXd x;
  double f = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
  
  EXPECT_FLOAT_EQ(6.4, f);
  EXPECT_NEAR(0, x(0), 1e-12);  // Need to deal with floating point precision
  EXPECT_FLOAT_EQ(2, x(1));
  EXPECT_FLOAT_EQ(0, x(2));
  
  // cout << "f: " << solve_quadprog(G, g0, CE, ce0, CI, ci0, x) << "\n"
  //      << "x: " << x << "\n";
}


struct fh {
  template <typename T0>
  inline Eigen::Matrix<T0, Eigen::Dynamic, Eigen::Dynamic> 
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
               const Eigen::VectorXd& delta) const {
      int n = 2;
      return Eigen::MatrixXd::Identity(n, n);
    }
};

struct fv {
  template <typename T0>
  inline Eigen::Matrix<T0, Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
               const Eigen::VectorXd& delta) const {
      int n = 2;
      Eigen::Matrix<T0, Eigen::Dynamic, 1> linear_term(n);
      linear_term(0) = 0;
      linear_term(1) = theta(1);
      return linear_term;
    }
};

struct fa_0 {
  template <typename T0>
  Eigen::Matrix<T0, Eigen::Dynamic, Eigen::Dynamic>
  inline operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                    const Eigen::VectorXd& delta) const {
    int n = 2;
    Eigen::Matrix<T0, Eigen::Dynamic, Eigen::Dynamic> linear_constraint(n, n);
    linear_constraint(0, 0) = 0;
    linear_constraint(0, 1) = 0;
    linear_constraint(1, 0) = 0;
    linear_constraint(1, 1) = 0;
    return linear_constraint;
  }
};

struct fa {
  template <typename T0>
  Eigen::Matrix<T0, Eigen::Dynamic, Eigen::Dynamic>
  inline operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                  const Eigen::VectorXd& delta) const {
    int n = 2;
    Eigen::Matrix<T0, Eigen::Dynamic, Eigen::Dynamic> linear_constraint(n, n);
    linear_constraint(0, 0) = 1;
    linear_constraint(0, 1) = 0;
    linear_constraint(1, 0) = 0;
    linear_constraint(1, 1) = 0;
    return linear_constraint;
  }
};

struct fb {
  template <typename T0>
  Eigen::Matrix<T0, Eigen::Dynamic, 1>
  inline operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                    const Eigen::VectorXd& delta) const {
    Eigen::Matrix<T0, Eigen::Dynamic, 1> b(2);
    b(0) = theta(0);
    b(1) = 0;
    return b;
  }
};

TEST(MathMatrix, quadratic_optimizer) {
  VectorXd theta(2);
  theta << 0, -1;
  VectorXd delta;
  int n = 2;

  VectorXd x = quadratic_optimizer(fh(), fv(), fa_0(), fb(), theta, delta, n);
  EXPECT_EQ(x(0), 0);
  EXPECT_EQ(x(1), 1);
  
  theta << -1, -1;
  x = quadratic_optimizer(fh(), fv(), fa(), fb(), theta, delta, n);
  EXPECT_EQ(x(0), 1);
  EXPECT_EQ(x(1), 1);
  
  // note: if constraint forces theta = -1, then inequality constraint
  // is overritten. Not sure if this should be concerning.
  // This can be seen by writing the test with theta << 1, -1,
  // which returns x = {-1, 1}.
  
  // test the inequality constraint.
  theta << 0, 1;
  x = quadratic_optimizer(fh(), fv(), fa_0(), fb(), theta, delta, n);
  std::cout << x << "\n \n";
}



