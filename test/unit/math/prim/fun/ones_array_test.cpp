#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, ones_array) {
  for (int K = 0; K < 5; K++) {
    std::vector<double> v(K, 1);
    EXPECT_STD_VECTOR_FLOAT_EQ(v, stan::math::ones_array(K));
  }
}

TEST(MathFunctions, ones_array_throw) {
  EXPECT_THROW(stan::math::ones_array(-1), std::domain_error);
}
