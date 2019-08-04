#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, lambert_w1) {
  using stan::math::lambert_w;
  using stan::math::is_nan;

  auto test_val = lambert_w(10.0, -1);
  printf("lambert_w(10): %.16f %s", test_val, "\n");
  EXPECT_TRUE(is_nan(test_val));

  test_val = lambert_w(-1.0, -1);
  printf("lambert_w(-1): %.16f %s", test_val, "\n");
  EXPECT_TRUE(is_nan(test_val));

  test_val = lambert_w(1.0, -1);
  printf("lambert_w(1): %.16f %s", test_val, "\n");
  EXPECT_TRUE(is_nan(test_val));

  test_val = lambert_w(-0.3, -1);
  printf("lambert_w(-0.3): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-1.7813370234216279, test_val);

  test_val = lambert_w(-0.2, -1);
  printf("lambert_w(-0.2): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-2.5426413577735261, test_val);
  // Test close to bound
  test_val = lambert_w(-0.367879, -1);
  printf("lambert_w(-0.367879): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-1.0015494951912602417, test_val);

  test_val = lambert_w(-0.367, -1);
  printf("lambert_w(-0.367): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-1.0707918867680526, test_val);

  test_val = lambert_w(-0.000001, -1);
  printf("lambert_w(-0.000001): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-16.62650890137247, test_val);

}

TEST(MathFunctions, lambert_w0) {
  using stan::math::lambert_w;
  using stan::math::is_nan;

  auto test_val = lambert_w(-1);
  printf("lambert_w(-1): %.16f %s", test_val, "\n");
  EXPECT_TRUE(is_nan(test_val));

  test_val = lambert_w(10.0);
  printf("lambert_w(10): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(1.7455280027406996, test_val);

  test_val = lambert_w(-0.3);
  printf("lambert_w(-0.3): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-0.4894022271802149, test_val);

  test_val = lambert_w(1.0);
  printf("lambert_w(1): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(0.5671432904097840, test_val);


  test_val = lambert_w(-0.2);
  printf("lambert_w(-0.2): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-0.2591711018190738, test_val);

  // Test close to bound
  test_val = lambert_w(-0.367879);
  printf("lambert_w(-0.367879): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-0.99845210378074245039, test_val);

  test_val = lambert_w(-0.367);
  printf("lambert_w(-0.367): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(-0.93239918474792793379, test_val);

  test_val = lambert_w(0.04);
  printf("lambert_w(0.04): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(0.038489665941978570829, test_val);
  // Test branch 0 for large values
  test_val = lambert_w(10000);
  printf("lambert_w(10000): %.16f %s", test_val, "\n");
  EXPECT_DOUBLE_EQ(7.231846038093372, test_val);

}
