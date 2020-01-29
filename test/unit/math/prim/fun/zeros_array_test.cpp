#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, zeros_array) {
  for (int K = 0; K < 5; K++) {
    std::vector<double> u = stan::math::zeros_array(K);
    std::vector<double> v(K, 0);
    EXPECT_EQ(u.size(), v.size());
    for (int i = 0; i < K; i++) {
      EXPECT_FLOAT_EQ(u[i], v[i]);
    }
  }
}

TEST(MathFunctions, zeros_array_throw) {
  EXPECT_THROW(stan::math::zeros_array(-1), std::domain_error);
}
