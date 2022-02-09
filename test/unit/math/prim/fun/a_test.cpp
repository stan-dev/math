#include <stan/math/prim.hpp>
#include <stan/math/prim/fun/a.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, a) {
  auto val_func = [&](auto&& x1, auto&& x2, auto&& x3) { return x1 + x2 + x3 + 1; };
  auto grad_func = [&](auto&& x1, auto&& x2, auto&& x3) { return x1 + x2 + x3 - 1; };
  stan::math::user_gradients(1, 2, 3, val_func, grad_func);
}
