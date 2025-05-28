#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

TEST(MathMetaPrim, is_tuple) {
  EXPECT_TRUE((stan::math::is_tuple<std::tuple<double>>::value));
  EXPECT_TRUE((stan::math::is_tuple<std::tuple<double>&>::value));
  EXPECT_TRUE((stan::math::is_tuple<std::tuple<double>&&>::value));
  EXPECT_FALSE((stan::math::is_tuple<double>::value));
  EXPECT_FALSE((stan::math::is_tuple<std::vector<double>>::value));
}
