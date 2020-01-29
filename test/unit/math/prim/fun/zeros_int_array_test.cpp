#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, zeros_int_array) {
  for (int K = 0; K < 5; K++) {
    std::vector<int> u = stan::math::zeros_int_array(K);
    std::vector<int> v(K, 0);
    EXPECT_EQ(u.size(), v.size());
    for (int i = 0; i < K; i++) {
      EXPECT_EQ(u[i], v[i]);
    }
  }
}

TEST(MathFunctions, zeros_int_array_throw) {
  EXPECT_THROW(stan::math::zeros_int_array(-1), std::domain_error);
}
