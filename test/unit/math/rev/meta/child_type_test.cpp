#include <stan/math/rev/meta.hpp>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>
#include <gtest/gtest.h>

TEST(MathMetaRevScal, value_type) {
  using stan::math::child_type;
  using stan::math::var;

  expect_same_type<double, child_type<var>::type>();
}
