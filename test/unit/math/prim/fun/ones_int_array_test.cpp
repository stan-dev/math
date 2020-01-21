#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, ones_int_array) {
  using stan::math::ones_int_array;

  std::vector<int> u0 = ones_int_array(0);
  EXPECT_EQ(0, u0.size());

  for (int K = 1; K < 5; K++) {
    std::vector<int> u = ones_int_array(K);
    std::vector<int> v(K, 1);
    EXPECT_EQ(u.size(), v.size());
    for (int i = 0; i < K; i++) {
      EXPECT_EQ(u[i], v[i]);
    }
  }
}

TEST(MathFunctions, ones_int_array_throw) {
  EXPECT_THROW(stan::math::ones_int_array(-1), std::domain_error);
}
