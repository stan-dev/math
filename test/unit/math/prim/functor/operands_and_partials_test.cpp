#include <stan/math/prim/functor.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, OperandsAndPartials) {
  using stan::math::operands_and_partials;

  auto o1 = operands_and_partials(1.0);
  auto o2 = operands_and_partials(2.0, 3.0, 4.0, 5.0);

  // TODO(Sean): Learn why it's size 6
  // EXPECT_GT(sizeof(&o2), sizeof(o2));

  EXPECT_FLOAT_EQ(27.1, o1.build(27.1));
}
