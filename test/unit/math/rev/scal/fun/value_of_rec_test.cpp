#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, value_of_rec) {
  using stan::math::var;
  using stan::math::value_of_rec;


  var v_a(5.0);

  EXPECT_FLOAT_EQ(5.0, value_of_rec(v_a));
}
