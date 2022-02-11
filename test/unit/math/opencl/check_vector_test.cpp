#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandlingOpenCL, checkVectorMatrix) {
  using stan::math::matrix_cl;
  matrix_cl<double> x(3, 3);

  EXPECT_THROW(stan::math::check_vector("checkVector", "x", x),
               std::invalid_argument);
  matrix_cl<double> x1(0, 0);
  EXPECT_THROW(stan::math::check_vector("checkVector", "x", x1),
               std::invalid_argument);

  matrix_cl<double> x2(1, 5);
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x2));

  matrix_cl<double> x3(5, 1);
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x3));
}

TEST(ErrorHandlingOpenCL, checkVectorMatrix_nan) {
  using stan::math::matrix_cl;
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> nan_vec({nan, nan, nan, nan, nan, nan, nan, nan, nan});
  matrix_cl<double> x(nan_vec, 3, 3);
  EXPECT_THROW(stan::math::check_vector("checkVector", "x", x),
               std::invalid_argument);

  matrix_cl<double> x1(nan_vec, 1, 9);
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x1));

  matrix_cl<double> x2(nan_vec, 9, 1);
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x2));
}
#endif
