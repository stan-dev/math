#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(ErrorHandlingMatrix, checkNonzeroSize) {
  using stan::math::check_nonzero_size;
  double y(0);
  
  EXPECT_TRUE(check_nonzero_size("checkNonzeroSize",
                                 "y", y));
}

TEST(ErrorHandlingMatrix, checkNonzeroSizeMatrix_nan) {
  using stan::math::check_nonzero_size;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(check_nonzero_size("checkNonzeroSize",
                                 "y", nan));
}
