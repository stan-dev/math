#include <stan/math/rev.hpp>
#include <stan/math/mix.hpp>
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
using stan::math::pow;

/*
 * A set of functions to test the finite-difference approximation on.
 */

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
    return stan::math::exp(2 * x(0)) + stan::math::exp(x(1));
  }
};

// exp_full(x, y) = e^{x-y}
struct exp_full {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return stan::math::exp(x(0) - x(1));
  }
};

// one_arg(x) = x^3
struct one_arg {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    using stan::math::pow;
    return pow(x(0), 3);
  }
};

/*
 * Comparing the finite-diff approximation to the AD hessian.
 */

template <typename F>
void test_hessian_finite_diff(const std::string& msg, const F& f,
                              Eigen::VectorXd& x) {
  double fx;
  Eigen::VectorXd grad_fx;
  Eigen::MatrixXd hess_fx;
  stan::math::internal::finite_diff_hessian_auto(f, x, fx, grad_fx, hess_fx);

  double fx_ad;
  Eigen::VectorXd grad_fx_ad;
  Eigen::MatrixXd hess_fx_ad;
  stan::math::hessian(f, x, fx_ad, grad_fx_ad, hess_fx_ad);

  EXPECT_FLOAT_EQ(fx_ad, fx) << msg;

  EXPECT_EQ(grad_fx_ad.size(), grad_fx.size());
  for (int i = 0; i < grad_fx_ad.size(); ++i)
    EXPECT_NEAR(grad_fx_ad(i), grad_fx(i), 1e-5) << msg;

  EXPECT_EQ(hess_fx_ad.rows(), hess_fx.rows()) << msg;
  EXPECT_EQ(hess_fx_ad.cols(), hess_fx.cols()) << msg;
  for (int i = 0; i < hess_fx_ad.size(); ++i)
    EXPECT_NEAR(hess_fx_ad(i), hess_fx(i), 1e-4) << msg;
}

TEST(RevFunctor, polynomial) {
  poly f;
  Matrix<double, Dynamic, 1> x(2);
  x << 5, 7;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_hessian_auto(f, x, fx, grad_fx, hess_fx);
  EXPECT_FLOAT_EQ(5 * 5 * 7 + 3 * 7 * 7, fx);
  EXPECT_EQ(2, grad_fx.size());
  EXPECT_EQ(2, hess_fx.rows());
  EXPECT_EQ(2, hess_fx.cols());
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
  stan::math::internal::finite_diff_hessian_auto(f, x, fx, grad_fx, hess_fx);
  EXPECT_EQ(3, hess_fx.rows());
  EXPECT_EQ(3, hess_fx.cols());
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
  stan::math::internal::finite_diff_hessian_auto(f, x, fx, grad_fx, hess_fx);
  EXPECT_EQ(2, hess_fx.rows());
  EXPECT_EQ(2, hess_fx.cols());
  EXPECT_FLOAT_EQ(4 * stan::math::exp(2 * 2), hess_fx(0, 0));
  EXPECT_FLOAT_EQ(0, hess_fx(0, 1));
  EXPECT_FLOAT_EQ(0, hess_fx(1, 0));
  EXPECT_FLOAT_EQ(stan::math::exp(-1), hess_fx(1, 1));
}

TEST(RevFunctor, exp_full) {
  exp_full f;
  Matrix<double, Dynamic, 1> x(2);
  x << 1, -3;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_hessian_auto(f, x, fx, grad_fx, hess_fx);
  EXPECT_EQ(2, hess_fx.rows());
  EXPECT_EQ(2, hess_fx.cols());
  EXPECT_FLOAT_EQ(stan::math::exp(4), hess_fx(0, 0));
  EXPECT_FLOAT_EQ(-stan::math::exp(4), hess_fx(0, 1));
  EXPECT_FLOAT_EQ(-stan::math::exp(4), hess_fx(1, 0));
  EXPECT_FLOAT_EQ(stan::math::exp(4), hess_fx(1, 1));
}

TEST(RevFunctor, one_arg) {
  one_arg f;
  Matrix<double, Dynamic, 1> x(1);
  x << 8;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  Matrix<double, Dynamic, Dynamic> hess_fx;
  stan::math::internal::finite_diff_hessian_auto(f, x, fx, grad_fx, hess_fx);
  EXPECT_EQ(1, hess_fx.rows());
  EXPECT_EQ(1, hess_fx.cols());
  EXPECT_FLOAT_EQ(6 * 8, hess_fx(0, 0));
}

TEST(RevFunctor, FiniteDiffHessianAuto) {
  auto norm_fun
      = [](const auto& x) { return stan::math::normal_lpdf(x(0), x(1), x(2)); };
  Eigen::VectorXd x(3);
  for (const auto& scale : std::vector<double>{1, 1e10, 1e20, 1e30}) {
    x << 1 * scale, 2 * scale, 1.5 * scale;
    test_hessian_finite_diff("norm_fun({1, 2, 3})", norm_fun, x);
  }
  for (const auto& scale : std::vector<double>{1e-10, 1e-20, 1e-30}) {
    x << scale, -scale, 1;  // finite diff fails with small sigma
    test_hessian_finite_diff("norm_fun({1, 2, 3})", norm_fun, x);
  }

  auto log_fun
      = [](const auto& x) { return stan::math::sum(stan::math::log(x)); };
  Eigen::VectorXd y(0);
  test_hessian_finite_diff("log_fun({})", log_fun, y);

  Eigen::VectorXd z(1);
  z << 2;
  test_hessian_finite_diff("log_fun({2})", log_fun, z);

  Eigen::VectorXd w(5);
  w << 1, 2, 3, 4, 5;
  test_hessian_finite_diff("log_fun({1, 2, 3, 4, 5})", log_fun, w);
}
