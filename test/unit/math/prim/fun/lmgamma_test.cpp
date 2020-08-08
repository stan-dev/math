#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, lmgamma) {
  unsigned int k = 1;
  double x = 2.5;
  double result = k * (k - 1) * log(boost::math::constants::pi<double>()) / 4.0;
  // j = 1
  result += lgamma(x);
  EXPECT_FLOAT_EQ(result, stan::math::lmgamma(k, x));

  k = 2;
  x = 3.0;
  result = k * (k - 1) * log(boost::math::constants::pi<double>()) / 4.0;
  // j = 1
  result += lgamma(x);
  // j = 2
  result += lgamma(x + (1.0 - 2.0) / 2.0);
  EXPECT_FLOAT_EQ(result, stan::math::lmgamma(k, x));
}

TEST(MathFunctions, lmgamma_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::lmgamma(2, nan)));
}

TEST(MathFunctions, lmgamma_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::lmgamma;
    return lmgamma(x1, x2);
  };

  std::vector<int> std_in1{1, 3, 1};
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 3);
  stan::test::binary_scalar_tester(f, std_std_in1, mat_in2);
}
