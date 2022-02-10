#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <stan/math/rev/fun/a.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, a) {
  stan::math::var a = 1;

  auto val_func = [&](auto&& x1) { return x1 + 1; };
  auto grad_func = [&](auto&& x1) { return x1 - 1; };
  stan::math::user_gradients(std::make_tuple(a),val_func,std::make_tuple(grad_func));
}
