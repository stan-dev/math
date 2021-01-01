#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandlingOpenCL, checkColumnIndexMatrix) {
  stan::math::matrix_cl<double> y;
  size_t i;

  i = 2;
  y = stan::math::matrix_cl<double>(3, 3);
  EXPECT_NO_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i));
  i = 3;
  EXPECT_NO_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i));

  y = stan::math::matrix_cl<double>(3, 2);
  EXPECT_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i),
      std::out_of_range);

  i = 0;
  EXPECT_THROW(
      stan::math::check_column_index("checkColumnIndexMatrix", "i", y, i),
      std::out_of_range);
}

#endif
