#include <stan/math/prim.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsDiscreteRange, error_check) {
  using stan::math::discrete_range_rng;
  boost::random::mt19937 rng;

  std::vector<int> lower{5, 11, -15};
  std::vector<int> upper{7, 15, -10};

  EXPECT_NO_THROW(discrete_range_rng(lower, upper, rng));
  EXPECT_THROW(discrete_range_rng(lower, 10, rng), std::domain_error);
  EXPECT_THROW(discrete_range_rng(10, upper, rng), std::domain_error);

  std::vector<int> vec2{1, 2};
  EXPECT_THROW(discrete_range_rng(lower, vec2, rng), std::invalid_argument);
  EXPECT_THROW(discrete_range_rng(vec2, upper, rng), std::invalid_argument);
}

TEST(ProbDistributionsDiscreteRange, boundary_values) {
  using stan::math::discrete_range_rng;
  boost::random::mt19937 rng;

  std::vector<int> lower{-5, 11, 17};
  EXPECT_EQ(lower, discrete_range_rng(lower, lower, rng));

  std::vector<int> upper(lower);
  for (int i = 0; i < upper.size(); i++) {
    ++upper[i];
  }

  EXPECT_LE(lower, discrete_range_rng(lower, upper, rng));
  EXPECT_GE(upper, discrete_range_rng(lower, upper, rng));
}

TEST(ProbDistributionsDiscreteRange, chiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;

  int N = 10000;
  int lower = -3;
  int upper = 20;

  int K = upper - lower + 1;
  boost::math::chi_squared mydist(K - 1);

  std::vector<int> bin(K, 0);
  std::vector<double> expect(K, static_cast<double>(N) / K);

  for (int count = 0; count < N; ++count) {
    int a = stan::math::discrete_range_rng(lower, upper, rng);
    ++bin[a - lower];
  }

  double chi = 0;

  for (int j = 0; j < K; j++) {
    chi += (bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j];
  }

  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}
