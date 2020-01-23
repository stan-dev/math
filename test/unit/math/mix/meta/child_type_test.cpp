#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/fun/promote_type_test_util.hpp>

TEST(MathMetaMix, value_type) {
  using stan::math::child_type;
  using stan::math::fvar;
  using stan::math::var;

  expect_same_type<fvar<var>, child_type<fvar<fvar<var> > >::type>();
  expect_same_type<var, child_type<fvar<var> >::type>();
}
