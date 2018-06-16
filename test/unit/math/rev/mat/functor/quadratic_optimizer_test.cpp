#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/quadratic_optimizer.hpp>
#include <stan/math/rev/mat/functor/eiquadprog.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

using stan::math::quadratic_optimizer;
using stan::math::quadratic_optimizer_analytical;
using stan::math::f_theta;
using stan::math::var;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;
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
  Eigen::Matrix<T0, Eigen::Dynamic, 1>
  inline operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                    const Eigen::VectorXd& delta) const {
    int n = 2;
    Eigen::Matrix<T0, Eigen::Dynamic, 1> linear_constraint(n);
    linear_constraint(0) = 0;
    linear_constraint(1) = 0;
    return linear_constraint;
  }
};

struct fa {
  template <typename T0>
  Eigen::Matrix<T0, Eigen::Dynamic, 1>
  inline operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                  const Eigen::VectorXd& delta) const {
    int n = 2;
    Eigen::Matrix<T0, Eigen::Dynamic, 1> linear_constraint(n);
    linear_constraint(0) = 1;
    linear_constraint(1) = 0;
    return linear_constraint;
  }
};

struct fb {
  template <typename T0>
  T0
  inline operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                    const Eigen::VectorXd& delta) const {
    return theta(0);
  }
};

// template <typename F1, typename F2, typename F3, typename F4>
// struct f_theta {
//   VectorXd delta_;
//   VectorXd x_;
// 
//   // structure to be passed to the stan::math::jacobian function
//   f_theta(const F1& fh,
//           const F2& fv,
//           const F3& fa,
//           const F4& fb,
//           const VectorXd& delta,
//           const VectorXd& x) : delta_(delta), x_(x) {}
// 
//   template<typename T>
//   Eigen::Matrix<T, Eigen::Dynamic, 1>
//   inline operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
//     return quadratic_optimizer_analytical(fh(), fv(), fa(), fb(),
//                                           theta, delta_, x_);
//   }
// };

TEST(MathMatrix, quadratic_optimizer_double) {
  VectorXd theta(2);
  theta << 0, -1;
  VectorXd delta;
  int n = 2;

  VectorXd x = quadratic_optimizer(fh(), fv(), fa_0(), fb(), theta, delta, n);
  EXPECT_EQ(x(0), 0);
  EXPECT_EQ(x(1), 1);

  theta << -5, -3;
  x = quadratic_optimizer(fh(), fv(), fa(), fb(), theta, delta, n);
  EXPECT_EQ(x(0), -theta(0));
  EXPECT_EQ(x(1), -theta(1));
  
  // Test analytical solution
  VectorXd x_an
    = quadratic_optimizer_analytical(fh(), fv(), fa(), fb(),
                                     theta, delta, x);
  EXPECT_EQ(x_an(0), -theta(0));
  EXPECT_EQ(x_an(1), -theta(1));

  MatrixXd J_solution(2, 2);
  J_solution << -1, 0,
                0, -1;
  
  // Test gradient computation with analytical solution
  f_theta<fh, fv, fa, fb> f(fh(), fv(), fa(), fb(), delta, x);
  VectorXd f_eval;
  MatrixXd J;
  stan::math::jacobian(f, theta, f_eval, J);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) EXPECT_EQ(J_solution(i, j), J(i, j));

    // cout << J << "\n";
  
  // theta << -1, 0;
  // std::cout << "Finite diff test \n"
  //           << quadratic_optimizer(fh(), fv(), fa(), fb(), theta, delta, n)
  //           << "\n \n";
  // 
  // double finite_diff = 1e-10;  // , finite_diff_2 = 1e-10;
  // VectorXd theta_lb = theta, theta_ub = theta;
  // theta_lb(1) = theta(1) - finite_diff;
  // theta_ub(1) = theta(1) + finite_diff;
  // 
  // VectorXd
  //   x_lb = quadratic_optimizer(fh(), fv(), fa(), fb(), theta_lb, delta, n),
  //   x_ub = quadratic_optimizer(fh(), fv(), fa(), fb(), theta_ub, delta, n);
  // 
  // std::cout << x_lb << "\n \n";
  // std::cout << x_ub << "\n \n";
  // 
  // std::cout << "dx_2 / d_theta_2 = "
  //           << (x_lb(1) - x_ub(1)) / (2 * finite_diff) << "\n \n";

  // note: if constraint forces theta = -1, then inequality constraint
  // is overritten. Not sure if this should be concerning.
  // This can be seen by writing the test with theta << 1, -1,
  // which returns x = {-1, 1}.

  // test the inequality constraint.
  // theta << 0, 1;
  // x = quadratic_optimizer(fh(), fv(), fa_0(), fb(), theta, delta, n);
  // std::cout << "Inequality constraint test:\n" << x << "\n \n";
}

TEST(MathMatrix, quadratic_optimizer_var) {
  VectorXd delta;
  int n_x = 2;
  int n_theta = 2;
  
  MatrixXd J(2, 2);
  J << -1, 0, 
       0, -1;

  for (int k = 0; k < n_x; k++) {
    Matrix<var, Dynamic, 1> theta(n_theta);
    theta << -5, -3;
    Matrix<var, Dynamic, 1> x = quadratic_optimizer(fh(), fv(), fa(), fb(),
                                                    theta, delta, n_x);
    EXPECT_EQ(-theta(0), x(0).val());
    EXPECT_EQ(-theta(1), x(1).val());
  
    AVEC parameter_vec = createAVEC(theta(0), theta(1));
    VEC g;
    x(k).grad(parameter_vec, g);
    EXPECT_EQ(J(k, 0), g[0]);
    EXPECT_EQ(J(k, 1), g[1]);
  }
}



