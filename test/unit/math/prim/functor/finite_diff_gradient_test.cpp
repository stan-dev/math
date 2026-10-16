#include <stan/math/prim/functor/finite_diff_gradient.hpp>
#include <gtest/gtest.h>

TEST(MathPrimFunctor, FiniteDiffGradientConstant) {
  // Preserve the evaluations when testing cancellation under FMA contraction.
  const volatile double constant = -937.6;
  auto f = [&constant](const Eigen::VectorXd&) { return constant; };
  Eigen::VectorXd x(3);
  x << 0, 1.7, -1000;
  double fx;
  Eigen::VectorXd grad;
  stan::math::finite_diff_gradient(f, x, fx, grad);

  EXPECT_DOUBLE_EQ(-937.6, fx);
  ASSERT_EQ(x.size(), grad.size());
  for (int i = 0; i < grad.size(); ++i) {
    EXPECT_DOUBLE_EQ(0, grad(i));
  }
}

TEST(MathPrimFunctor, FiniteDiffGradientPolynomial) {
  auto f = [](const Eigen::VectorXd& x) {
    return x(0) * x(0) * x(0) + x(1) * x(1) * x(1) * x(1);
  };
  Eigen::VectorXd x(2);
  x << 1.7, -0.3;
  double fx;
  Eigen::VectorXd grad;
  stan::math::finite_diff_gradient(f, x, fx, grad);

  EXPECT_DOUBLE_EQ(f(x), fx);
  ASSERT_EQ(2, grad.size());
  EXPECT_NEAR(3 * x(0) * x(0), grad(0), 1e-10);
  EXPECT_NEAR(4 * x(1) * x(1) * x(1), grad(1), 1e-10);
}
