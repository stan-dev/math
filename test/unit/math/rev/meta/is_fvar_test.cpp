#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, is_fvar) {
  using stan::is_fvar;
  EXPECT_FALSE(is_fvar<stan::math::var>::value);
}
