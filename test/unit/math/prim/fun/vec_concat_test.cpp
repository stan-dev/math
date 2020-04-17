#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, vec_concat) {
  std::vector<int> x{0, 1, 2, 3};
  auto x_result = stan::math::vec_concat(x, std::vector<int>{4, 5, 6},
                                         std::vector<int>{7, 8, 9});
  EXPECT_TRUE((x_result.size() == 10));
  for (int i = 0; i < 10; i++) {
    EXPECT_EQ(i, x_result[i]);
  }
}
