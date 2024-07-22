#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, checkStochasticColumn) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_vec(2, 2);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y{
      y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i << 0.5, 0.5, 0.5, 0.5;
  }

  EXPECT_NO_THROW(
      stan::math::check_stochastic_column("checkStochasticColumn", "y", y));

  for (auto& y_i : y) {
    y_i(0, 1) = 0.55;
  }
  EXPECT_THROW(
      stan::math::check_stochastic_column("checkStochasticColumn", "y", y),
      std::domain_error);
}

TEST(ErrorHandlingMatrix, checkStochasticColumn_message_negative_value) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_vec(3, 3);
  y_vec.setZero();
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y{
      y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i(0, 0) = -0.1;
    y_i(1, 0) = 1.1;
    y_i(0, 1) = -0.1;
    y_i(1, 1) = 1.1;
  }

  try {
    stan::math::check_stochastic_column("checkStochasticColumn", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_TRUE(std::string::npos
              != message.find(" y[1] is not a valid column stochastic matrix"))
      << "Found: " << message;

  EXPECT_TRUE(std::string::npos != message.find("y[1][1, 1] = -0.1"))
      << "Found: " << message;

  for (auto& y_i : y) {
    y_i.setZero();
    y_i(0, 0) = 0.1;
    y_i(1, 0) = 0.1;
    y_i(2, 0) = 1.0;
  }
  try {
    stan::math::check_stochastic_column("checkStochasticColumn", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_TRUE(std::string::npos
              != message.find(" y[1] is not a valid column stochastic matrix"))
      << "Found: " << message;

  EXPECT_TRUE(std::string::npos != message.find("sum(y[1][:, 1]) = 1.2"))
      << "Found: " << message;
}

TEST(ErrorHandlingMatrix, checkStochasticColumn_message_sum) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_vec(10, 10);
  y_vec.setZero();
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y{
      y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i(3, 0) = 0.9;
  }

  try {
    stan::math::check_stochastic_column("checkStochasticColumn", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_TRUE(std::string::npos
              != message.find(" y[1] is not a valid column stochastic matrix"))
      << message;

  EXPECT_TRUE(std::string::npos != message.find("sum(y[1][:, 1]) = 0.9"))
      << message;
}

TEST(ErrorHandlingMatrix, checkStochasticColumn_message_length) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_vec(0, 0);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y{
      y_vec, y_vec, y_vec};

  using stan::math::check_stochastic_column;

  EXPECT_THROW_MSG(check_stochastic_column("checkStochasticColumn", "y", y),
                   std::invalid_argument,
                   "y[1] has size 0, but must have a non-zero size");
}

TEST(ErrorHandlingMatrix, checkStochasticColumn_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_vec(2, 2);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y{
      y_vec, y_vec, y_vec};
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  for (auto& y_i : y) {
    y_i << nan, 0.5, nan, 0.5;
  }

  EXPECT_THROW(
      stan::math::check_stochastic_column("checkStochasticColumn", "y", y),
      std::domain_error);

  for (auto& y_i : y) {
    y_i(0, 1) = 0.55;
  }
  EXPECT_THROW(
      stan::math::check_stochastic_column("checkStochasticColumn", "y", y),
      std::domain_error);

  for (auto& y_i : y) {
    y_i(0, 0) = 0.5;
    y_i(0, 1) = nan;
  }
  EXPECT_THROW(
      stan::math::check_stochastic_column("checkStochasticColumn", "y", y),
      std::domain_error);

  for (auto& y_i : y) {
    y_i(0, 0) = nan;
  }
  EXPECT_THROW(
      stan::math::check_stochastic_column("checkStochasticColumn", "y", y),
      std::domain_error);
}
