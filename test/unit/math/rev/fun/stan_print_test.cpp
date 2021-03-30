#include <stan/math/rev.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, stan_print) {
  using stan::math::var;

  var a = 5.0;

  {
    std::stringstream s;
    stan::math::stan_print(&s, a);
    EXPECT_TRUE(s.str().find("5") != std::string::npos);
  }
}
