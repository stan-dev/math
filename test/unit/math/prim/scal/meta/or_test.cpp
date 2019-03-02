#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/meta/and.hpp>
#include <gtest/gtest.h>

TEST(MathMeta, or_type) {
    bool temp = stan::math::or_<std::true_type, std::true_type, std::true_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_<std::false_type, std::false_type, std::false_type>::value;
    EXPECT_FALSE(temp);
    temp = stan::math::or_<std::false_type, std::true_type, std::true_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_<std::true_type, std::true_type, std::false_type>::value;
    EXPECT_TRUE(temp);
}

TEST(MathMeta, or_not_type) {
    bool temp = stan::math::or_not_<std::false_type, std::false_type, std::false_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_not_<std::false_type, std::false_type, std::true_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_not_<std::false_type, std::true_type, std::false_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_not_<std::false_type, std::true_type, std::true_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_not_<std::true_type, std::false_type, std::false_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_not_<std::true_type, std::false_type, std::true_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_not_<std::true_type, std::true_type, std::false_type>::value;
    EXPECT_TRUE(temp);
    temp = stan::math::or_not_<std::true_type, std::true_type, std::true_type>::value;
    EXPECT_FALSE(temp);
}