#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <stan/math/prim/fun/a.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, a) {
  stan::math::var a = 1;
  stan::math::user_gradients(a, 2.5, "b", "c", "d");
}
