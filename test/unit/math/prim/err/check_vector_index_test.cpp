#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/err/util.hpp>
#include <limits>

TEST(ErrorHandlingMatrix, checkVectorIndexColumnVector) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  size_t i;
  double nan = std::numeric_limits<double>::quiet_NaN();

  i = 2;
  y.resize(3);
  y << nan, nan, nan;
  EXPECT_NO_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i));
  i = 3;
  EXPECT_NO_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i));

  y.resize(2);
  y << nan, nan;
  STAN_EXPECT_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i),
      std::out_of_range);

  i = 0;
  STAN_EXPECT_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i),
      std::out_of_range);
}

TEST(ErrorHandlingMatrix, checkVectorIndexRowVector) {
  Eigen::Matrix<double, 1, Eigen::Dynamic> y;
  size_t i;
  double nan = std::numeric_limits<double>::quiet_NaN();

  i = 2;
  y.resize(3);
  y << nan, nan, nan;
  EXPECT_NO_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i));
  i = 3;
  EXPECT_NO_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i));

  y.resize(2);
  y << nan, nan;
  STAN_EXPECT_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i),
      std::out_of_range);

  i = 0;
  STAN_EXPECT_THROW(
      stan::math::check_vector_index("checkVectorIndexMatrix", "i", y, i),
      std::out_of_range);
}
