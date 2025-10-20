#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, MetaTraitsRevScal_is_fvar) {
  using stan::is_fvar;
  EXPECT_FALSE(is_fvar<stan::math::var>::value);
}
