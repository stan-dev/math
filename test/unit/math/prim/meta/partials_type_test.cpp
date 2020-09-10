#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, PartialsTypeDouble) {
  test::expect_same_type<double, stan::partials_type<double>::type>();
}

TEST(MathMetaPrim, PartialsTypeFloat) {
  test::expect_same_type<float, stan::partials_type<float>::type>();
}

TEST(MathMetaPrim, PartialsTypeInt) {
  test::expect_same_type<int, stan::partials_type<int>::type>();
}
