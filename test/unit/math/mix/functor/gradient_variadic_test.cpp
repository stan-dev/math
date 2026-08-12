#include <stan/math/mix/functor/gradient_variadic.hpp>
#include <gtest/gtest.h>

namespace gradient_variadic_test {

struct cubic_functor {
  template <typename T0, typename T1, typename T2>
  auto operator()(const T0& theta, std::ostream* /* msgs */, const T1& b,
                  const T2& c) const {
    auto theta_cubed = theta(0) * theta(0) * theta(0);
    for (Eigen::Index i = 1; i < theta.size(); ++i) {
      theta_cubed += theta(i) * theta(i) * theta(i);
    }
    return c * theta_cubed / 6.0 + b.dot(theta);
  }
};

}  // namespace gradient_variadic_test

TEST(MathMixFunctor, GradientVariadicArithmetic) {
  Eigen::VectorXd theta(3);
  Eigen::VectorXd b(3);
  theta << 0.4, -0.7, 1.2;
  b << 0.3, -0.2, 0.5;
  const double c = 1.3;

  Eigen::VectorXd g = stan::math::gradient(
      gradient_variadic_test::cubic_functor{}, theta, nullptr, b, c);
  Eigen::VectorXd expected
      = 0.5 * c * theta.array().square().matrix() + b;
  EXPECT_TRUE(g.isApprox(expected, 1e-12));
}

TEST(MathMixFunctor, GradientVariadicReverse) {
  using stan::math::var;
  Eigen::VectorXd theta(3);
  Eigen::VectorXd b(3);
  Eigen::VectorXd w(3);
  theta << 0.4, -0.7, 1.2;
  b << 0.3, -0.2, 0.5;
  w << -0.5, 0.8, 0.2;
  const double c = 1.3;

  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var = theta;
  Eigen::Matrix<var, Eigen::Dynamic, 1> b_var = b;
  var c_var = c;
  auto g = stan::math::gradient(
      gradient_variadic_test::cubic_functor{}, theta_var, nullptr, b_var,
      c_var);
  var target = w.dot(g);
  target.grad();

  Eigen::VectorXd expected
      = 0.5 * c * theta.array().square().matrix() + b;
  Eigen::VectorXd theta_adj = c * w.array() * theta.array();
  EXPECT_TRUE(g.val().isApprox(expected, 1e-12));
  EXPECT_TRUE(theta_var.adj().isApprox(theta_adj, 1e-12));
  EXPECT_TRUE(b_var.adj().isApprox(w, 1e-12));
  EXPECT_NEAR(c_var.adj(),
              0.5 * (w.array() * theta.array().square()).sum(), 1e-12);
  stan::math::recover_memory();
}

TEST(MathMixFunctor, GradientVariadicVarArgument) {
  using stan::math::var;
  Eigen::VectorXd theta(3);
  Eigen::VectorXd b(3);
  Eigen::VectorXd w(3);
  theta << 0.4, -0.7, 1.2;
  b << 0.3, -0.2, 0.5;
  w << -0.5, 0.8, 0.2;
  var c = 1.3;

  auto g = stan::math::gradient(
      gradient_variadic_test::cubic_functor{}, theta, nullptr, b, c);
  w.dot(g).grad();
  EXPECT_NEAR(c.adj(), 0.5 * (w.array() * theta.array().square()).sum(),
              1e-12);
  stan::math::recover_memory();
}

TEST(MathMixFunctor, GradientVariadicEmpty) {
  Eigen::VectorXd theta(0);
  Eigen::VectorXd b(0);
  EXPECT_EQ(stan::math::gradient(
                gradient_variadic_test::cubic_functor{}, theta, nullptr, b,
                1.0)
                .size(),
            0);
}
