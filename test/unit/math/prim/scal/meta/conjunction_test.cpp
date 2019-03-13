#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <gtest/gtest.h>

TEST(MathMeta, and_type) {
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
