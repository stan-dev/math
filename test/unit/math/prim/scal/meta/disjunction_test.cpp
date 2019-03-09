#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <gtest/gtest.h>

TEST(MathMeta, or_type) {
  bool temp = stan::math::disjunction<std::true_type, std::true_type,
                                      std::true_type>::value;
  EXPECT_TRUE(temp);
  temp = stan::math::disjunction<std::false_type, std::false_type,
                                 std::false_type>::value;
  EXPECT_FALSE(temp);
  temp = stan::math::disjunction<std::false_type, std::true_type,
                                 std::true_type>::value;
  EXPECT_TRUE(temp);
  temp = stan::math::disjunction<std::true_type, std::true_type,
                                 std::false_type>::value;
  EXPECT_TRUE(temp);
}
