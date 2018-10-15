#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>

TEST(AgradRevMatrix, softmax) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::softmax;
  using stan::math::vector_v;

  EXPECT_THROW(softmax(vector_v()), std::invalid_argument);

  Matrix<AVAR, Dynamic, 1> x(1);
  x << 0.0;

  Matrix<AVAR, Dynamic, 1> theta = softmax(x);
  EXPECT_EQ(1, theta.size());
  EXPECT_FLOAT_EQ(1.0, theta[0].val());

  Matrix<AVAR, Dynamic, 1> x2(2);
  x2 << -1.0, 1.0;
  Matrix<AVAR, Dynamic, 1> theta2 = softmax(x2);
  EXPECT_EQ(2, theta2.size());
  EXPECT_FLOAT_EQ(exp(-1) / (exp(-1) + exp(1)), theta2[0].val());
  EXPECT_FLOAT_EQ(exp(1) / (exp(-1) + exp(1)), theta2[1].val());

  Matrix<AVAR, Dynamic, 1> x3(3);
  x3 << -1.0, 1.0, 10.0;
  Matrix<AVAR, Dynamic, 1> theta3 = softmax(x3);
  EXPECT_EQ(3, theta3.size());
  EXPECT_FLOAT_EQ(exp(-1) / (exp(-1) + exp(1) + exp(10.0)), theta3[0].val());
  EXPECT_FLOAT_EQ(exp(1) / (exp(-1) + exp(1) + exp(10.0)), theta3[1].val());
  EXPECT_FLOAT_EQ(exp(10) / (exp(-1) + exp(1) + exp(10.0)), theta3[2].val());
}

TEST(AgradRevSoftmax, gradient_check) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::softmax;
  using stan::math::var;
  std::vector<std::vector<double> > inputs
      = {{0.5, -1.0, 3.0}, {4.0, 3.0, -2.0}};
  std::vector<double> vals = {0.07459555713221443, 0.729736214118415};
  std::vector<std::vector<double> > grads
      = {{0.06903105998834897, -0.001241607138856182, -0.06778945284949277},
         {0.1972212719225377, -0.1959012993504623, -0.001319972572075391}};
  for (size_t i = 0; i < inputs.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      stan::math::set_zero_all_adjoints();

      Matrix<AVAR, Dynamic, 1> alpha(3);
      for (int k = 0; k < 3; ++k)
        alpha((j + k) % 3) = inputs[i][k];

      Matrix<AVAR, Dynamic, 1> theta = softmax(alpha);
      theta(j).grad();

      EXPECT_NEAR(vals[i], theta(j).val(), 1e-10);
      for (size_t k = 0; k < 3; ++k)
        EXPECT_NEAR(grads[i][k], alpha((j + k) % 3).adj(), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> alpha(3);
  alpha << 0.0, 3.0, -1.0;
  test::check_varis_on_stack(stan::math::softmax(alpha));
}
