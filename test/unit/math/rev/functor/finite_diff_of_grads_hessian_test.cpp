#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <stdexcept>
#include <vector>
#include <thread>
#include <future>

using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::exp;
using std::pow;

// poly(x, y) = (x^2 * y) + (3 * y^2)
struct poly {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1);
  }
};

// linear(x, y, z) = 3x-1y+8z
struct linear {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return 3 * x(0) - x(1) + 8 * x(2);
  }
};

// exp_diag(x, y) = e^2x+e^y
struct exp_diag {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return exp(2 * x(0)) + exp(x(1));
  }
};

// exp_full(x, y) = e^{x-y}
struct exp_full {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return exp(x(0) - x(1));
  }
};

// one_arg(x) = x^3
struct one_arg {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return pow(x(0), 3);
  }
};

TEST(RevFunctor, polynomial) {
  poly f;
  Matrix<double, Dynamic, 1> x(2);
  x << 5, 7;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_of_grads_hessian(f, x, fx, grad_fx,
                                                     hess_fx);
  EXPECT_FLOAT_EQ(5 * 5 * 7 + 3 * 7 * 7, fx);
  EXPECT_EQ(2, grad_fx.size());
  EXPECT_FLOAT_EQ(2 * x(0) * x(1), grad_fx(0));
  EXPECT_FLOAT_EQ(x(0) * x(0) + 3 * 2 * x(1), grad_fx(1));
  EXPECT_FLOAT_EQ(2 * x(1), hess_fx(0, 0));
  EXPECT_FLOAT_EQ(2 * x(0), hess_fx(0, 1));
  EXPECT_FLOAT_EQ(2 * x(0), hess_fx(1, 0));
  EXPECT_FLOAT_EQ(6, hess_fx(1, 1));
}

TEST(RevFunctor, linear_function) {
  linear f;
  Matrix<double, Dynamic, 1> x(3);
  x << 5, 7, -1;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_of_grads_hessian(f, x, fx, grad_fx,
                                                     hess_fx);
  EXPECT_FLOAT_EQ(0, hess_fx(0, 0));
  EXPECT_FLOAT_EQ(0, hess_fx(0, 1));
  EXPECT_FLOAT_EQ(0, hess_fx(0, 2));
  EXPECT_FLOAT_EQ(0, hess_fx(1, 0));
  EXPECT_FLOAT_EQ(0, hess_fx(1, 1));
  EXPECT_FLOAT_EQ(0, hess_fx(1, 2));
  EXPECT_FLOAT_EQ(0, hess_fx(2, 0));
  EXPECT_FLOAT_EQ(0, hess_fx(2, 1));
  EXPECT_FLOAT_EQ(0, hess_fx(2, 2));
}

TEST(RevFunctor, exp_diag) {
  exp_diag f;
  Matrix<double, Dynamic, 1> x(2);
  x << 2, -1;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_of_grads_hessian(f, x, fx, grad_fx,
                                                     hess_fx);
  EXPECT_FLOAT_EQ(4 * exp(2 * 2), hess_fx(0, 0));
  EXPECT_FLOAT_EQ(0, hess_fx(0, 1));
  EXPECT_FLOAT_EQ(0, hess_fx(1, 0));
  EXPECT_FLOAT_EQ(exp(-1), hess_fx(1, 1));
}

TEST(RevFunctor, exp_full) {
  exp_full f;
  Matrix<double, Dynamic, 1> x(2);
  x << 1, -3;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_of_grads_hessian(f, x, fx, grad_fx,
                                                     hess_fx);
  EXPECT_FLOAT_EQ(exp(4), hess_fx(0, 0));
  EXPECT_FLOAT_EQ(-exp(4), hess_fx(0, 1));
  EXPECT_FLOAT_EQ(-exp(4), hess_fx(1, 0));
  EXPECT_FLOAT_EQ(exp(4), hess_fx(1, 1));
}

TEST(RevFunctor, one_arg) {
  one_arg f;
  Matrix<double, Dynamic, 1> x(1);
  x << 8;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_of_grads_hessian(f, x, fx, grad_fx,
                                                     hess_fx);
  EXPECT_FLOAT_EQ(6 * 8, hess_fx(0, 0));
}