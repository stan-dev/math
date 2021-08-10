#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, PartialsTypeDouble) {
  EXPECT_SAME_TYPE(double, stan::partials_type<double>::type);
}

TEST(MathMetaPrim, PartialsTypeFloat) {
  EXPECT_SAME_TYPE(float, stan::partials_type<float>::type);
}

TEST(MathMetaPrim, PartialsTypeInt) {
  EXPECT_SAME_TYPE(int, stan::partials_type<int>::type);
}
