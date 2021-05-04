#include <stan/math/prim/meta.hpp>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, value_type) {
  using stan::math::child_type;

  expect_same_type<double, child_type<double>::type>();
}
