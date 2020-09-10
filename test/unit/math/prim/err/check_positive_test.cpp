#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandlingArr, CheckPositive) {
  using stan::math::check_positive;
  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    EXPECT_NO_THROW(check_positive("check_positive_test", "x", x));
  }
}

TEST(ErrorHandlingArr, CheckPositive_nan) {
  using stan::math::check_positive;
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = nan;
    EXPECT_THROW(check_positive("check_positive_test", "x", x),
                 std::domain_error);
    x[i] = i;
  }
}

TEST(ErrorHandlingMat, CheckPositive) {
  using stan::math::check_positive;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x_mat(3);
  x_mat << 1, 2, 3;
  for (int i = 0; i < x_mat.size(); i++) {
    EXPECT_NO_THROW(check_positive("check_positive_test", "x", x_mat));
  }

  x_mat(0) = 0;
  EXPECT_THROW(check_positive("check_positive_test", "x", x_mat),
               std::domain_error);
}

TEST(ErrorHandlingMat, CheckPositive_nan) {
  using stan::math::check_positive;
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double, Eigen::Dynamic, 1> x_mat(3);
  x_mat << 1, 2, 3;
  for (int i = 0; i < x_mat.size(); i++) {
    x_mat(i) = nan;
    EXPECT_THROW(check_positive("check_positive_test", "x", x_mat),
                 std::domain_error);
    x_mat(i) = i;
  }
}

TEST(ErrorHandlingScalar, CheckPositive) {
  using stan::math::check_positive;
  EXPECT_NO_THROW(check_positive("check_positive_test", "x", 3.0));
}

TEST(ErrorHandlingScalar, CheckPositive_nan) {
  using stan::math::check_positive;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_positive("check_positive_test", "x", nan),
               std::domain_error);
}

TEST(ErrorHandlingScalar, CheckPositive_0) {
  using stan::math::check_positive;
  EXPECT_THROW(check_positive("check_positive_test", "x", 0u),
               std::domain_error);
  EXPECT_THROW(check_positive("check_positive_test", "x", (size_t)0),
               std::domain_error);
  EXPECT_THROW(check_positive("check_positive_test", "x", 0.0),
               std::domain_error);
  EXPECT_THROW(check_positive("check_positive_test", "x", -0.0),
               std::domain_error);
  EXPECT_THROW(check_positive("check_positive_test", "x", 0),
               std::domain_error);
}

TEST(ErrorHandlingScalar, CheckPositiveVectorization) {
  using stan::math::check_positive;
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 1);
  EXPECT_NO_THROW(check_positive("check_positive_test", "m",
                                 std::vector<Eigen::MatrixXd>{m, m, m}));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = -1;
  EXPECT_THROW(check_positive("check_positive_test", "m",
                              std::vector<Eigen::MatrixXd>{m, m2, m}),
               std::domain_error);
}
