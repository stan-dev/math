#include <stan/math/prim.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsDiscreteRange, error_check) {
  using stan::math::discrete_range_rng;
  boost::random::mt19937 rng;

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  Eigen::VectorXd lower(3);
  Eigen::VectorXd upper(3);

  lower << 5.0, 11.0, -15.0;
  upper << 5.0, 15.0, -10.0;
  EXPECT_NO_THROW(discrete_range_rng(lower, upper, rng));
  EXPECT_THROW(discrete_range_rng(lower, 10, rng), std::domain_error);
  EXPECT_THROW(discrete_range_rng(10, upper, rng), std::domain_error);

  lower << -1e3, 1.1e3, 1e4;
  upper << -1e2, 1.2e3, 1e5;
  EXPECT_NO_THROW(discrete_range_rng(lower, upper, rng));

  EXPECT_THROW(discrete_range_rng(nan, upper, rng), std::domain_error);
  EXPECT_THROW(discrete_range_rng(inf, upper, rng), std::domain_error);
  EXPECT_THROW(discrete_range_rng(-inf, upper, rng), std::domain_error);
  EXPECT_THROW(discrete_range_rng(lower, nan, rng), std::domain_error);
  EXPECT_THROW(discrete_range_rng(lower, inf, rng), std::domain_error);
  EXPECT_THROW(discrete_range_rng(lower, -inf, rng), std::domain_error);

  Eigen::VectorXd vec2(2);
  vec2 << 1, 2;
  EXPECT_THROW(discrete_range_rng(lower, vec2, rng), std::invalid_argument);
  EXPECT_THROW(discrete_range_rng(vec2, upper, rng), std::invalid_argument);
}

TEST(ProbDistributionsDiscreteRange, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;

  int N = 10000;
  int lower = -3;
  int upper = 20;

  int K = upper - lower + 1;
  boost::math::chi_squared mydist(K - 1);

  int bin[K];
  double expect[K];
  double prop = static_cast<double>(N) / K;
  for (int i = 0; i < K; i++) {
    bin[i] = 0;
    expect[i] = prop;
  }

  int count = 0;
  while (count < N) {
    int a = stan::math::discrete_range_rng(lower, upper, rng);
    bin[a - lower]++;
    count++;
  }

  double chi = 0;

  for (int j = 0; j < K; j++) {
    chi += (bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j];
  }

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
