#include <gtest/gtest.h>
#include <sstream>
#include <stan/math/mix/arr.hpp>

TEST(AgradFwdFvar, insertion_operator) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a(5.0);
  std::stringstream ss;
  ss << a;
  EXPECT_EQ("5", ss.str());
}
