#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

const char* function = "function";
class ErrorHandlingMatrix : public ::testing::Test {
public:
  void SetUp() {
  }

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
};

TEST_F(ErrorHandlingMatrix, checkPosSemidefinite_size_1) {
  using stan::math::check_pos_semidefinite;

  y.resize(1, 1);

  y << 0.0;
  EXPECT_NO_THROW(check_pos_semidefinite(function, "y", y));

  y << -1.0;
  EXPECT_THROW_MSG(check_pos_semidefinite(function, "y", y),
                   std::domain_error,
                   "function: y is not positive semi-definite.");
}

TEST_F(ErrorHandlingMatrix, checkPosSemidefinite_bad_sizes) {
  using stan::math::check_pos_semidefinite;

  y.resize(0,0);
  EXPECT_THROW_MSG(check_pos_semidefinite(function, "y", y),
                   std::invalid_argument,
                   "function: y must have a positive size, but is 0; dimension size expression = rows");

  y.resize(2,3);
  EXPECT_THROW_MSG(check_pos_semidefinite(function, "y", y),
                   std::invalid_argument,
                   "function: Expecting a square matrix; rows of y (2) and columns of y (3) must match in size");
}

TEST_F(ErrorHandlingMatrix, checkPosSemidefinite) {
  using stan::math::check_pos_semidefinite;

  y.resize(2, 2);

  y << 1, 0, 0, 1;
  EXPECT_NO_THROW(check_pos_semidefinite(function, "y", y));

  y << -1, 0, 0, 1;
  EXPECT_THROW_MSG(check_pos_semidefinite(function, "y", y),
                   std::domain_error,
                   "function: y is not positive semi-definite.");
}

TEST_F(ErrorHandlingMatrix, checkPosSemidefinite_nan) {
  using stan::math::check_pos_semidefinite;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3,3);
  y << 2, -1, 0,
    -1, 2, -1,
    0, -1, 2;
  EXPECT_NO_THROW(check_pos_semidefinite(function, "y", y));
  for (int i = 0; i < y.rows(); i++)
    for (int j = 0; j < y.cols(); j++) {
      y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
      y(i,j) = nan;
      if (i >= j) {
        // std::stringstream expected_msg;
        // if (i == j) {
        //   expected_msg << "function: y[" << j * y.cols() + i + 1 << "] is "
        //                << nan << ", but must not be nan!";
        // } else {
        //   expected_msg << "function: y is not symmetric. "
        //                << "y[" << j + 1 << "," << i + 1 << "] = " << y(j, i)
        //                << ", but y[" << i + 1 << "," << j + 1 << "] = " << y(i, j);
        // }
        EXPECT_THROW(check_pos_semidefinite(function, "y", y),
                     std::domain_error);  // , expected_msg.str());
      }
    }

}
