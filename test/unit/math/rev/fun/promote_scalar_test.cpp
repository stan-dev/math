#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

// there is no agrad-defined version of promote_scalar, so this is
// just testing that it works with non-inter-convertible types (double
// can be assigned to var, but not vice-versa)

TEST(AgradRevFunctionsPromoteScalar, Mismatch) {
  using stan::math::promote_scalar;
  using stan::math::var;
  EXPECT_FLOAT_EQ(2.3, promote_scalar<var>(2.3).val());
}
