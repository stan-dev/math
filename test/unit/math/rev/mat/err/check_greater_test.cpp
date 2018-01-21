#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

using stan::math::check_greater;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar, CheckGreaterMatrix) {
  const char* function = "check_greater";
  var x;
  var low;
  Eigen::Matrix<var, Eigen::Dynamic, 1> x_vec;
  Eigen::Matrix<var, Eigen::Dynamic, 1> low_vec;
  x_vec.resize(3);
  low_vec.resize(3);

  // x_vec, low_vec
  x_vec << -1, 0, 1;
  low_vec << -2, -1, 0;
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>";

  x_vec << -1, 0, 1;
  low_vec << -1.1, -0.1, 0.9;
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>";

  x_vec << -1, 0, std::numeric_limits<double>::infinity();
  low_vec << -2, -1, 0;
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>, y has infinity";

  x_vec << -1, 0, 1;
  low_vec << -2, 0, 0;
  EXPECT_THROW(check_greater(function, "x", x_vec, low_vec), std::domain_error)
      << "check_greater: matrix<3, 1>, matrix<3, 1>, should pass for index 1";

  x_vec << -1, 0, 1;
  low_vec << -2, -1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x_vec, low_vec), std::domain_error)
      << "check_greater: matrix<3, 1>, matrix<3, 1>, should fail with infinity";

  x_vec << -1, 0, std::numeric_limits<double>::infinity();
  low_vec << -2, -1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x_vec, low_vec), std::domain_error)
      << "check_greater: matrix<3, 1>, matrix<3, 1>, "
      << "both bound and value infinity";

  x_vec << -1, 0, 1;
  low_vec << -2, -1, -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>, should pass with "
         "-infinity";

  // x_vec, low
  x_vec << -1, 0, 1;
  low = -2;
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low))
      << "check_greater: matrix<3, 1>, double";

  x_vec << -1, 0, 1;
  low = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low))
      << "check_greater: matrix<3, 1>, double";

  x_vec << -1, 0, 1;
  low = 0;
  EXPECT_THROW(check_greater(function, "x", x_vec, low), std::domain_error)
      << "check_greater: matrix<3, 1>, double, should fail for index 1/2";

  x_vec << -1, 0, 1;
  low = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x_vec, low), std::domain_error)
      << "check_greater: matrix<3, 1>, double, should fail with infinity";

  // x, low_vec
  x = 2;
  low_vec << -1, 0, 1;
  EXPECT_NO_THROW(check_greater(function, "x", x, low_vec))
      << "check_greater: double, matrix<3, 1>";

  x = 0;
  low_vec << -1, 0, -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x, low_vec), std::domain_error)
      << "check_greater: double, matrix<3, 1>, low has -inf";

  x = 10;
  low_vec << -1, 0, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x, low_vec), std::domain_error)
      << "check_greater: double, matrix<3, 1>, low has inf";

  x = std::numeric_limits<double>::infinity();
  low_vec << -1, 0, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x, low_vec), std::domain_error)
      << "check_greater: double, matrix<3, 1>, x is inf, low has inf";

  x = std::numeric_limits<double>::infinity();
  low_vec << -1, 0, 1;
  EXPECT_NO_THROW(check_greater(function, "x", x, low_vec))
      << "check_greater: double, matrix<3, 1>, x is inf";

  x = 1.1;
  low_vec << -1, 0, 1;
  EXPECT_NO_THROW(check_greater(function, "x", x, low_vec))
      << "check_greater: double, matrix<3, 1>";

  x = 0.9;
  low_vec << -1, 0, 1;
  EXPECT_THROW(check_greater(function, "x", x, low_vec), std::domain_error)
      << "check_greater: double, matrix<3, 1>";
  stan::math::recover_memory();
}
