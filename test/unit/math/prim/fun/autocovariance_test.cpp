#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/ar1.hpp>
#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <stdexcept>

TEST(ProbAutocovariance, test1) {
  // ar1.csv generated in R with
  //   > x[1] <- rnorm(1, 0, 1)
  //   > for (n in 2:1000) x[n] <- rnorm(1, 0.8 * x[n-1], 1)
   std::vector<double> y{stan::math::test::ar1_arr.begin(),
    stan::math::test::ar1_arr.end()};

  // 10K 1K-length AC in 2.9s with g++ -O3 on Bob's Macbook Air
  std::vector<double> ac;
  stan::math::autocovariance(y, ac);

  EXPECT_EQ(1000U, ac.size());

  EXPECT_NEAR(2.69, ac[0], 0.01);
  EXPECT_NEAR(2.16, ac[1], 0.01);
  EXPECT_NEAR(1.73, ac[2], 0.01);
  EXPECT_NEAR(1.38, ac[3], 0.01);
  EXPECT_NEAR(1.10, ac[4], 0.01);
  EXPECT_NEAR(0.90, ac[5], 0.01);
}

TEST(ProbAutocovariance, test2) {
  static constexpr size_t N = 1000;
  Eigen::VectorXd y(N);
  for (size_t i = 0; i < N; ++i) {
    y(i) = stan::math::test::ar1_arr[i];
  }

  // 10K 1K-length AC in 2.9s with g++ -O3 on Bob's Macbook Air
  Eigen::VectorXd ac(N);
  stan::math::autocovariance<double>(y, ac);

  EXPECT_EQ(1000U, ac.size());

  EXPECT_NEAR(2.69, ac(0), 0.01);
  EXPECT_NEAR(2.16, ac(1), 0.01);
  EXPECT_NEAR(1.73, ac(2), 0.01);
  EXPECT_NEAR(1.38, ac(3), 0.01);
  EXPECT_NEAR(1.10, ac(4), 0.01);
  EXPECT_NEAR(0.90, ac(5), 0.01);
}
