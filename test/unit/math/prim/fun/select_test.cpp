#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, select_scalar_scalar) {
  using stan::math::select;

  EXPECT_FLOAT_EQ(2.0, select(true, 2.0, 5.5));
  EXPECT_FLOAT_EQ(5.5, select(false, 2.0, 5.5));

  double dbl_a = -5.5;
  int int_a = -5;
  double dbl_b = 6.3251;
  int int_b = 6;

  int int_res = select(false, dbl_a, int_b);
  double dbl_res = select(true, dbl_a, int_b);

  EXPECT_EQ(int_res, int_b);
  EXPECT_FLOAT_EQ(dbl_res, dbl_a);
}

TEST(MathFunctions, select_container_container) {
  using stan::math::select;

  std::vector<double> std_res{1.2, 63, 1.25};

  EXPECT_STD_VECTOR_EQ(2.0, select(true, 2.0, 5.5));
  EXPECT_FLOAT_EQ(5.5, select(false, 2.0, 5.5));

  double dbl_a = -5.5;
  int int_a = -5;
  double dbl_b = 6.3251;
  int int_b = 6;

  int int_res = select(false, dbl_a, int_b);
  double dbl_res = select(true, dbl_a, int_b);

  EXPECT_EQ(int_res, int_b);
  EXPECT_FLOAT_EQ(dbl_res, dbl_a);
}
