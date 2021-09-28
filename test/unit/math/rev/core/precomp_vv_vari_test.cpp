#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>

TEST(StanAgradRevInternal, precomp_vv_vari) {
  double value, gradient1, gradient2;
  stan::math::var x1(2), x2(3);
  stan::math::var y;

  value = 1;
  gradient1 = 4;
  gradient2 = 5;

  std::vector<stan::math::var> vars{x1, x2};

  EXPECT_NO_THROW(y = stan::math::var(new stan::math::precomp_vv_vari(
                      value, x1.vi_, x2.vi_, gradient1, gradient2)));
  EXPECT_FLOAT_EQ(value, y.val());

  std::vector<double> g;
  EXPECT_NO_THROW(y.grad(vars, g));
  ASSERT_EQ(2U, g.size());
  EXPECT_FLOAT_EQ(gradient1, g[0]);
  EXPECT_FLOAT_EQ(gradient2, g[1]);

  stan::math::recover_memory();
}
