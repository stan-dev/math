#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, and_type) {
  bool temp = stan::math::conjunction<std::true_type, std::true_type,
                                      std::true_type>::value;
  EXPECT_TRUE(temp);
  temp = stan::math::conjunction<std::false_type, std::false_type,
                                 std::false_type>::value;
  EXPECT_FALSE(temp);
  temp = stan::math::conjunction<std::false_type, std::true_type,
                                 std::true_type>::value;
  EXPECT_FALSE(temp);
  temp = stan::math::conjunction<std::true_type, std::true_type,
                                 std::false_type>::value;
  EXPECT_FALSE(temp);
}
