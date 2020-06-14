#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <sstream>

TEST(AgradRev, insertion_operator) {
  stan::math::vari v(5);
  std::stringstream ss;
  ss << &v;
  EXPECT_EQ("5:0", ss.str());
}

TEST(AgradRev, long_double_test) {
  stan::math::vari_value<float> v(5);
  std::stringstream ss;
  ss << &v;
  EXPECT_EQ("5:0", ss.str());
}
