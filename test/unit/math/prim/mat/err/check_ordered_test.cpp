#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using stan::math::check_ordered;

TEST(ErrorHandlingMatrix, checkOrdered) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  y.resize(3);

  y << 0, 1, 2;
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y));

  y << 0, 10, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y));

  y << -10, 10, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y));

  y << -std::numeric_limits<double>::infinity(), 10, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y));

  y << 0, 0, 0;
  EXPECT_THROW(check_ordered("check_ordered", "y", y),
               std::domain_error);

  y << 0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_ordered("check_ordered", "y", y),
               std::domain_error);


  y << -1, 3, 2;
  EXPECT_THROW(check_ordered("check_ordered", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkOrdered_one_indexed_message) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  y.resize(3);
  
  y << 0, 5, 1;
  try {
    check_ordered("check_ordered", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("element at 3"))
    << message;
}

TEST(ErrorHandlingMatrix, checkOrdered_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  y.resize(3);

  y << 0, 1, 2;
  for (int i = 0; i < y.size(); i++) {
    y[i] = nan;
    EXPECT_THROW(check_ordered("check_ordered", "y", y),
                 std::domain_error);
    y[i] = i;
  }
  for (int i = 0; i < y.size(); i++) {
    y << 0, 10, std::numeric_limits<double>::infinity();
    y[i] = nan;
    EXPECT_THROW(check_ordered("check_ordered", "y", y),
                 std::domain_error);
  }
}
