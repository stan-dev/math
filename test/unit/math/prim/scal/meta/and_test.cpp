#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/meta/and.hpp>
#include <gtest/gtest.h>

TEST(MathMeta, and_type) {
    bool temp = stan::math::and_<std::true_type, std::true_type, std::true_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::and_<std::false_type, std::false_type, std::false_type>::value;
    EXPECT_FALSE(temp);
    temp = stan::math::and_<std::false_type, std::true_type, std::true_type>::value;
    EXPECT_FALSE(temp);
    temp = stan::math::and_<std::true_type, std::true_type, std::false_type>::value;
    EXPECT_FALSE(temp);
}
