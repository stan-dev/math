#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkColumnIndexMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  size_t i;

  i = 2;
  y.resize(3, 3);
  EXPECT_NO_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i));
  i = 3;
  EXPECT_NO_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i));

  y.resize(3, 2);
  EXPECT_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i),
      std::out_of_range);

  i = 0;
  EXPECT_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i),
      std::out_of_range);
}

TEST(ErrorHandlingMatrix, checkColumnIndexMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  size_t i;

  i = 2;
  y.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i));
  i = 3;
  EXPECT_NO_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i));

  y.resize(3, 2);
  y << nan, nan, nan, nan, nan, nan;
  EXPECT_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i),
      std::out_of_range);

  i = 0;
  EXPECT_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i),
      std::out_of_range);
}
