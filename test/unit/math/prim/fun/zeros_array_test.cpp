#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, zeros_array) {
  for (int K = 0; K < 5; K++) {
    std::vector<double> v(K, 0);
    expect_std_vector_eq(v, stan::math::zeros_array(K));
  }
}

TEST(MathFunctions, zeros_array_throw) {
  EXPECT_THROW(stan::math::zeros_array(-1), std::domain_error);
}
