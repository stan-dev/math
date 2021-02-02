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

// fun1(x, y) = (x^2 * y) + (3 * y^2)
struct fun1 {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1);
  }
};

TEST(RevFunctor, finite_diff_hessian) {
  fun1 f;
  Matrix<double, Dynamic, 1> x(2);
  x << 5, 7;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::finite_diff_of_grads_hessian(f, x, fx, grad_fx, hess_fx);
  EXPECT_FLOAT_EQ(5 * 5 * 7 + 3 * 7 * 7, fx);
  EXPECT_EQ(2, grad_fx.size());
  EXPECT_FLOAT_EQ(2 * x(0) * x(1), grad_fx(0));
  EXPECT_FLOAT_EQ(x(0) * x(0) + 3 * 2 * x(1), grad_fx(1));
  EXPECT_FLOAT_EQ(2 * x(1), hess_fx(0,0));
  EXPECT_FLOAT_EQ(2 * x(0), hess_fx(0,1));
  EXPECT_FLOAT_EQ(2 * x(0), hess_fx(1,0));
  EXPECT_FLOAT_EQ(6, hess_fx(1,1));
}

