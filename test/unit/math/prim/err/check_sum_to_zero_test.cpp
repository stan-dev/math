#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, checkSumToZero_edges) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> zero(0);

  EXPECT_THROW(stan::math::check_sum_to_zero("checkSumToZero", "zero", zero),
               std::invalid_argument);

  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(1);
  y_vec << 0.0;

  EXPECT_NO_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y_vec));

  y_vec[0] = 0.1;
  EXPECT_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y_vec),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkSumToZero) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(2);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i << 0.5, -0.5;
  }

  EXPECT_NO_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y));

  for (auto& y_i : y) {
    y_i[0] = 0.55;
  }
  EXPECT_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkSumToZero_message_sum) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(100);
  y_vec.setZero();
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i[13] = 0.1;
  }

  try {
    stan::math::check_sum_to_zero("checkSumToZero", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_TRUE(std::string::npos != message.find(" y[1] does not sum to zero"))
      << message;

  EXPECT_TRUE(std::string::npos != message.find("sum(y[1]) = 0.1")) << message;
}

TEST(ErrorHandlingMatrix, checkSumToZero_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(2);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  for (auto& y_i : y) {
    y_i << nan, 0.5;
  }

  EXPECT_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i[1] = 0.55;
  }
  EXPECT_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i[0] = 0.5;
    y_i[1] = nan;
  }
  EXPECT_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i[0] = nan;
  }
  EXPECT_THROW(stan::math::check_sum_to_zero("checkSumToZero", "y", y),
               std::domain_error);
}
