#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

using stan::math::check_ordered;

TEST(ErrorHandling, checkOrdered) {
  std::vector<double> y_;
  y_.push_back(0.0);
  y_.push_back(1.0);
  y_.push_back(2.0);
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y_));

  y_[2] = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y_));

  y_[0] = -10.0;
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y_));

  y_[0] = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_ordered("check_ordered", "y", y_));

  y_[0] = 0.0;
  y_[1] = 0.0;
  y_[2] = 0.0;
  EXPECT_THROW(check_ordered("check_ordered", "y", y_),
               std::domain_error);

  y_[1] = std::numeric_limits<double>::infinity();
  y_[2] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_ordered("check_ordered", "y", y_),
               std::domain_error);
}

TEST(ErrorHandling, checkOrdered_one_indexed_message) {
  std::string message;
  std::vector<double> y;
  y.push_back(0);
  y.push_back(5);
  y.push_back(1);
  
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

TEST(ErrorHandling, checkOrdered_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> y_;
  y_.push_back(0.0);
  y_.push_back(1.0);
  y_.push_back(2.0);
  for (size_t i = 0; i < y_.size(); i++) {
    y_[i] = nan;
    EXPECT_THROW(check_ordered("check_ordered", "y", y_),
                 std::domain_error);
    y_[i] = i;
  }
  for (size_t i = 0; i < y_.size(); i++) {
    y_[0] = 0.0;
    y_[1] = 10.0;
    y_[2] = std::numeric_limits<double>::infinity();
    y_[i] = nan;
    EXPECT_THROW(check_ordered("check_ordered", "y", y_),
                 std::domain_error);
  }
}
