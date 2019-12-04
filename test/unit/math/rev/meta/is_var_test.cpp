#include <stan/math/rev/meta.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, is_var) {
  using stan::is_var;
  EXPECT_TRUE(is_var<stan::math::var>::value);
}
