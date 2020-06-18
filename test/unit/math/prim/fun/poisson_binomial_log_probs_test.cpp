#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctionsPoissonBinomialLogProbs, check_scalar_y) {
  using stan::math::poisson_binomial_log_probs;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<double> vals{-2.12026, -0.84397, -0.967584, -2.65926};

  for (int i = 0; i < 4; i++) {
    auto x = poisson_binomial_log_probs(i, p);
    EXPECT_EQ(x.rows(), i + 1);
    EXPECT_EQ(x.cols(), 1);
    for (int j = 0; j <= i; ++j) {
      EXPECT_NEAR(vals[j], x(j), 0.001);
    }
  }
}

TEST(MathFunctionsPoissonBinomialLogProbs, check_vectorial_y) {
  using stan::math::poisson_binomial_log_probs;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{0, 1};
  auto x = poisson_binomial_log_probs(y, p);

  EXPECT_EQ(x.size(), 2);
  EXPECT_EQ(x[0].rows(), 1);
  EXPECT_EQ(x[0].cols(), 1);
  EXPECT_EQ(x[1].rows(), 2);
  EXPECT_EQ(x[1].cols(), 1);
}

TEST(MathFunctionsPoissonBinomialLogProbs, check_vectorial_y_vectorial_theta) {
  using stan::math::poisson_binomial_log_probs;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{0, 1};
  std::vector<vec> ps{p, p};
  auto x = poisson_binomial_log_probs(y, ps);

  EXPECT_EQ(x.size(), 2);
  EXPECT_EQ(x[0].rows(), 1);
  EXPECT_EQ(x[0].cols(), 1);
  EXPECT_EQ(x[1].rows(), 2);
  EXPECT_EQ(x[1].cols(), 1);
}
