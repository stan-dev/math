#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using stan::math::check_equal;

TEST(ErrorHandlingScalar,CheckEqualMatrix) {
  const char* function = "check_equal";
  Eigen::Matrix<double,Eigen::Dynamic,1> x_vec;
  Eigen::Matrix<double,Eigen::Dynamic,1> eq_vec;
  x_vec.resize(3);
  eq_vec.resize(3);

  // x_vec, low_vec
  x_vec   << -1, 0, 1;
  eq_vec << -1, 0, 1;
  EXPECT_TRUE(check_equal(function, "x", x_vec, eq_vec)) 
    << "check_equal: matrix<3,1>, matrix<3,1>";

  x_vec   <<   -1,    0,   1;
  eq_vec << -1.1, -0.1, 0.9;
  EXPECT_THROW(check_equal(function, "x", x_vec, eq_vec),
               std::domain_error) 
    << "check_equal: matrix<3,1>, matrix<3,1>";
  
  x_vec   << -1, 0,  1;
  eq_vec << -2, -1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_equal(function, "x", x_vec, eq_vec), 
               std::domain_error) 
    << "check_equal: matrix<3,1>, matrix<3,1>, should fail with infinity";
}


TEST(ErrorHandlingScalar,CheckEqual_Matrix_one_indexed_message) {
  const char* function = "check_equal";
  double x;
  double eq;
  Eigen::Matrix<double,Eigen::Dynamic,1> x_vec(3);
  Eigen::Matrix<double,Eigen::Dynamic,1> eq_vec(3);
  std::string message;
  
  // x_vec, eq_vec
  x_vec   <<   0,    0,   1;
  eq_vec <<    0,    0,   0;

  try {
    check_equal(function, "x", x_vec, eq_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }
  
  EXPECT_NE(std::string::npos, message.find("[3]"))
    << message;

  // x, eq_vec
  x = 1;
  try {
    check_equal(function, "x", x, eq_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }
  
  EXPECT_EQ(std::string::npos, message.find("["))
    << "no index information" << std::endl
    << message;

  // x_vec, eq
  eq = 1;
  try {
    check_equal(function, "x", x_vec, eq);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }
  
  EXPECT_NE(std::string::npos, message.find("[1]"))
    << message;

}

TEST(ErrorHandlingScalar,CheckEqual_nan) {
  const char* function = "check_equal";
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double,Eigen::Dynamic,1> x_vec;
  Eigen::Matrix<double,Eigen::Dynamic,1> eq_vec;
  x_vec.resize(3);
  eq_vec.resize(3);

  // x_vec, low_vec
  x_vec   << nan, 0, 1;
  eq_vec << -1, 0, 1;
  EXPECT_THROW(check_equal(function, "x", x_vec, eq_vec),
               std::domain_error);

  eq_vec << nan, 0, 1;
  EXPECT_THROW(check_equal(function, "x", x_vec, eq_vec),
               std::domain_error);

  x_vec << -1, 0, 1;
  EXPECT_THROW(check_equal(function, "x", x_vec, eq_vec),
               std::domain_error);
}
