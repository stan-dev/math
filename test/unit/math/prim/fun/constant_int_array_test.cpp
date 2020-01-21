#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, constant_int_array) {
  using stan::math::constant_int_array;

  std::vector<int> u0 = constant_int_array(0, 1);
  EXPECT_EQ(0, u0.size());

  int c = 8;
  for (int K = 1; K < 5; K++) {
    std::vector<int> u = constant_int_array(K, c);
    std::vector<int> v(K, c);
    EXPECT_EQ(u.size(), v.size());
    for (int i = 0; i < K; i++) {
      EXPECT_EQ(u[i], v[i]);
    }
  }
}

TEST(MathFunctions, constant_int_array_throw) {
  EXPECT_THROW(stan::math::constant_int_array(-1, 3), std::domain_error);
}
