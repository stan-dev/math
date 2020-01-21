#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MathFunctions, spaced_array) {
  using stan::math::spaced_array;

  std::vector<double> u0 = spaced_array(0, 1, 2);
  EXPECT_EQ(0, u0.size());

  double low = 1;
  double high = 5;
  std::vector<double> u11 = spaced_array(1, low, high);
  EXPECT_EQ(1, u11.size());
  EXPECT_FLOAT_EQ(high, u11[0]);

  int K = 5;
  std::vector<double> u1 = spaced_array(K, 1, 5);
  std::vector<double> v1{1, 2, 3, 4, 5};

  EXPECT_EQ(K, u1.size());
  for (int i = 0; i < K; i++) {
    EXPECT_FLOAT_EQ(u1[i], v1[i]);
  }

  std::vector<double> u2 = spaced_array(K, -2, 2);
  std::vector<double> v2{-2, -1, 0, 1, 2};

  EXPECT_EQ(K, u2.size());
  for (int i = 0; i < K; i++) {
    EXPECT_FLOAT_EQ(u2[i], v2[i]);
  }
}

TEST(MathFunctions, spaced_array_throw) {
  using stan::math::spaced_array;
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  int K = 5;
  double low = -2;
  double high = 6;

  EXPECT_THROW(spaced_array(-1, low, high), std::domain_error);

  EXPECT_THROW(spaced_array(K, inf, high), std::domain_error);
  EXPECT_THROW(spaced_array(K, nan, high), std::domain_error);

  EXPECT_THROW(spaced_array(K, low, low - 1), std::domain_error);
  EXPECT_THROW(spaced_array(K, low, inf), std::domain_error);
  EXPECT_THROW(spaced_array(K, low, nan), std::domain_error);
}
