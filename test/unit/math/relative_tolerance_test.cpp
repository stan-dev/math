#include <test/unit/math/relative_tolerance.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(testUnitMath, RelativeTolerance) {
  using stan::test::relative_tolerance;

  // Note that default tol_min is 1e-14 in this case
  relative_tolerance t_8(1e-8);
  EXPECT_FLOAT_EQ(t_8.exact(1e3), 1e-5);
  EXPECT_FLOAT_EQ(t_8.exact(1), 1e-8);
  EXPECT_FLOAT_EQ(t_8.exact(1e-6), 1e-14);
  EXPECT_FLOAT_EQ(t_8.exact(1e-14), 1e-14);
  EXPECT_FLOAT_EQ(t_8.exact(0), 1e-14);

  EXPECT_FLOAT_EQ(t_8.inexact(-1e3, 3e3), 2e-5);
  EXPECT_FLOAT_EQ(t_8.inexact(1, 1000), 5.005e-6);
  EXPECT_FLOAT_EQ(t_8.inexact(1e-5, -8e-6), 9e-14);
  EXPECT_FLOAT_EQ(t_8.inexact(0, 2e3), 1e-5);
  EXPECT_FLOAT_EQ(t_8.inexact(0, 1e-10), 1e-14);

  // Default tol_min is 1e-8
  relative_tolerance t_4(1e-4);
  EXPECT_FLOAT_EQ(t_4.exact(1e3), 1e-1);
  EXPECT_FLOAT_EQ(t_4.exact(1), 1e-4);
  EXPECT_FLOAT_EQ(t_4.exact(-1e-3), 1e-7);
  EXPECT_FLOAT_EQ(t_4.exact(1e-5), 1e-8);
  EXPECT_FLOAT_EQ(t_4.exact(0), 1e-8);

  EXPECT_FLOAT_EQ(t_4.inexact(-1e3, 3e3), 2e-1);
  EXPECT_FLOAT_EQ(t_4.inexact(2e-5, 0), 1e-8);
  EXPECT_FLOAT_EQ(t_4.inexact(0, 1), 5e-5);
  EXPECT_FLOAT_EQ(t_4.inexact(0, -1e-10), 1e-8);

  relative_tolerance t_16_10(1e-16, 1e-10);
  EXPECT_FLOAT_EQ(t_16_10.exact(-30), 1e-10);
  EXPECT_FLOAT_EQ(t_16_10.exact(-500), 1e-10);
  EXPECT_FLOAT_EQ(t_16_10.exact(1e-3), 1e-10);
  EXPECT_FLOAT_EQ(t_16_10.exact(-1e-5), 1e-10);

  EXPECT_FLOAT_EQ(t_16_10.inexact(-5e6, 1e-14), 2.5e-10);
  EXPECT_FLOAT_EQ(t_16_10.inexact(1e-315, 1e-2), 1e-10);
}
