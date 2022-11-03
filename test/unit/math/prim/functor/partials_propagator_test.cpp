#include <stan/math/prim/functor.hpp>
#include <gtest/gtest.h>

TEST(MathMetaPrim, PartialsPropagator) {
  using stan::math::make_partials_propagator;

  auto o1 = stan::math::make_partials_propagator(1.0);
  auto o2 = stan::math::make_partials_propagator(2.0, 3.0, 4.0, 5.0);

  // This is size 8 because it just hold an tuple with 4 elements of size 8
  EXPECT_EQ(8, sizeof(o2));

  EXPECT_FLOAT_EQ(27.1, o1.build(27.1));
}
