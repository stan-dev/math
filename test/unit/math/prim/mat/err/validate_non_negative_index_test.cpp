#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathPrimMat, val_non_neg_index) {
  const char* var_name = "var_name";
  const char* expr_name = "expr";
  const int val_good = 1;
  const int val_bad = -1;

  EXPECT_THROW(
      stan::math::validate_non_negative_index(var_name, expr_name, val_bad),
      std::invalid_argument);
  EXPECT_NO_THROW(
      stan::math::validate_non_negative_index(var_name, expr_name, val_good));
}
