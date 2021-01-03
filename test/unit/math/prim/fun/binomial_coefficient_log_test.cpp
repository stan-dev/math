#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <cmath>

template <typename T_N, typename T_n>
void test_binom_coefficient(const T_N& N, const T_n& n) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(lgamma(N + 1) - lgamma(n + 1) - lgamma(N - n + 1),
                  binomial_coefficient_log(N, n))
      << "N = " << N << ", n = " << n;
}

TEST(MathFunctions, binomial_coefficient_log) {
  using stan::math::binomial_coefficient_log;
  EXPECT_FLOAT_EQ(1.0, exp(binomial_coefficient_log(2.0, 2.0)));
  EXPECT_FLOAT_EQ(2.0, exp(binomial_coefficient_log(2.0, 1.0)));
  EXPECT_FLOAT_EQ(3.0, exp(binomial_coefficient_log(3.0, 1.0)));
  EXPECT_NEAR(3.0, exp(binomial_coefficient_log(3.0, 2.0)), 0.0001);

  EXPECT_FLOAT_EQ(29979.16, binomial_coefficient_log(100000, 91116));

  EXPECT_EQ(binomial_coefficient_log(-1, 0), 0);  // Needed for neg_binomial_2
  EXPECT_EQ(binomial_coefficient_log(50, 0), 0);
  EXPECT_EQ(binomial_coefficient_log(10000, 0), 0);

  EXPECT_EQ(binomial_coefficient_log(10, 11), stan::math::NEGATIVE_INFTY);
  EXPECT_EQ(binomial_coefficient_log(10, -1), stan::math::NEGATIVE_INFTY);

  for (int n = 0; n < 1010; ++n) {
    test_binom_coefficient(1010, n);
    test_binom_coefficient(1010.0, n);
    test_binom_coefficient(1010, static_cast<double>(n));
    test_binom_coefficient(1010.0, static_cast<double>(n));
  }

  test_binom_coefficient(1e9, 1e5);
  test_binom_coefficient(1e50, 1e45);
  test_binom_coefficient(1e20, 1e15);
}

TEST(MathFunctions, binomial_coefficient_log_nan) {
  double nan = stan::math::NOT_A_NUMBER;

  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(2.0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(nan, 2.0)));
  EXPECT_TRUE(std::isnan(stan::math::binomial_coefficient_log(nan, nan)));
}

TEST(MathFunctions, binomial_coefficient_log_errors_edge_cases) {
  using stan::math::binomial_coefficient_log;
  using stan::math::INFTY;

  EXPECT_NO_THROW(binomial_coefficient_log(10, 11));
  EXPECT_THROW(binomial_coefficient_log(10, 11.01), std::domain_error);
  EXPECT_THROW(binomial_coefficient_log(10, -1.1), std::domain_error);
  EXPECT_THROW(binomial_coefficient_log(-1, 0.3), std::domain_error);
  EXPECT_NO_THROW(binomial_coefficient_log(-0.5, 0.49));
  EXPECT_NO_THROW(binomial_coefficient_log(10, -0.9));

  EXPECT_FLOAT_EQ(binomial_coefficient_log(0, -1), -INFTY);
  EXPECT_FLOAT_EQ(binomial_coefficient_log(-1, 0), 0);
  EXPECT_FLOAT_EQ(binomial_coefficient_log(-1, -0.3), INFTY);
  EXPECT_FLOAT_EQ(binomial_coefficient_log(0.3, -1), -INFTY);
  EXPECT_FLOAT_EQ(binomial_coefficient_log(5.0, 6.0), -INFTY);
}

TEST(MathFunctions, binomial_coefficient_log_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::binomial_coefficient_log;
    return binomial_coefficient_log(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 6.5, 13, 15;
  Eigen::VectorXd in2(3);
  in2 << 0.2, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
