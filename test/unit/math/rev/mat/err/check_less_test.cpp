#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

using stan::math::check_less;
using stan::math::var;

TEST(AgradRevErrorHandlingScalar,CheckLess_Matrix) {
  const char* function = "check_less";
  var x;
  var high;
  Eigen::Matrix<var,Eigen::Dynamic,1> x_vec;
  Eigen::Matrix<var,Eigen::Dynamic,1> high_vec;
  x_vec.resize(3);
  high_vec.resize(3);
  
  
  // x_vec, high
  x_vec << -5, 0, 5;
  high = 10;
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high));

  x_vec << -5, 0, 5;
  high = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high));

  x_vec << -5, 0, 5;
  high = 5;
  EXPECT_THROW(check_less(function, "x", x_vec, high),
               std::domain_error);
  
  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high = 5;
  EXPECT_THROW(check_less(function, "x", x_vec, high),
               std::domain_error);

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_less(function, "x", x_vec, high),
               std::domain_error);
  
  // x_vec, high_vec
  x_vec << -5, 0, 5;
  high_vec << 0, 5, 10;
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high_vec));

  x_vec << -5, 0, 5;
  high_vec << std::numeric_limits<double>::infinity(), 10, 10;
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high_vec));

  x_vec << -5, 0, 5;
  high_vec << 10, 10, 5;
  EXPECT_THROW(check_less(function, "x", x_vec, high_vec),
               std::domain_error);
  
  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high_vec << 10, 10, 10;
  EXPECT_THROW(check_less(function, "x", x_vec, high_vec), 
               std::domain_error);

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high_vec << 10, 10, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_less(function, "x", x_vec, high_vec),
               std::domain_error);

  
  // x, high_vec
  x = -100;
  high_vec << 0, 5, 10;
  EXPECT_NO_THROW(check_less(function, "x", x, high_vec));

  x = 10;
  high_vec << 100, 200, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less(function, "x", x, high_vec));

  x = 5;
  high_vec << 100, 200, 5;
  EXPECT_THROW(check_less(function, "x", x, high_vec),
               std::domain_error);
  
  x = std::numeric_limits<double>::infinity();
  high_vec << 10, 20, 30;
  EXPECT_THROW(check_less(function, "x", x, high_vec), 
               std::domain_error);

  x = std::numeric_limits<double>::infinity();
  high_vec << std::numeric_limits<double>::infinity(), 
    std::numeric_limits<double>::infinity(), 
    std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_less(function, "x", x, high_vec),
               std::domain_error);
  stan::math::recover_memory();
}
