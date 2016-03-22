#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>

TEST(MathMeta, value_type) {
  using stan::math::child_type;
  using stan::math::var;

  expect_same_type<double,child_type<var>::type>();
}
