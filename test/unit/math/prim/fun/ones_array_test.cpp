#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, ones_array) {
  using stan::math::ones_array;

  std::vector<double> u0 = ones_array(0);
  EXPECT_EQ(0, u0.size());

  for (int K = 1; K < 5; K++) {
    std::vector<double> u = ones_array(K);
    std::vector<double> v(K, 1);
    EXPECT_EQ(u.size(), v.size());
    for (int i = 0; i < K; i++) {
      EXPECT_FLOAT_EQ(u[i], v[i]);
    }
  }
}

TEST(MathFunctions, ones_array_throw) {
  EXPECT_THROW(stan::math::ones_array(-1), std::domain_error);
}
