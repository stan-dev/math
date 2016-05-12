#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/scal/fun/promote_type_test_util.hpp>

TEST(MathMeta, value_type) {
  using stan::math::child_type;
  using stan::math::fvar;

  expect_same_type<double,child_type<fvar<double> >::type>();
  expect_same_type<fvar<double>,child_type<fvar<fvar<double> > >::type>();
}
