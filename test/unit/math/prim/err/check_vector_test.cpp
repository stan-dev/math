#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkVectorMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;

  x.resize(3, 3);
  EXPECT_THROW(stan::math::check_vector("checkVector", "x", x),
               std::invalid_argument);
  x.resize(0, 0);
  EXPECT_THROW(stan::math::check_vector("checkVector", "x", x),
               std::invalid_argument);

  x.resize(1, 5);
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x));

  x.resize(5, 1);
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x));
}

TEST(ErrorHandlingMatrix, checkVectorMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x;
  double nan = std::numeric_limits<double>::quiet_NaN();

  x.resize(3, 3);
  x << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_THROW(stan::math::check_vector("checkVector", "x", x),
               std::invalid_argument);
  x.resize(0, 0);
  EXPECT_THROW(stan::math::check_vector("checkVector", "x", x),
               std::invalid_argument);

  x.resize(1, 5);
  x << nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x));

  x.resize(5, 1);
  x << nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_vector("checkVector", "x", x));
}
