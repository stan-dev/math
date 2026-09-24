#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, zip_index) {
  const int J = 5;  // Number of locations
  const int K = 4;  // Number of time categories

  Eigen::MatrixXd delta(J, K);
  delta << 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44, 51,
      52, 53, 54;

  std::vector<int> time = {1, 1, 1, 2, 3, 3, 4, 1, 1};

  std::vector<int> loc = {1, 1, 5, 2, 3, 3, 4, 4, 5};

  Eigen::VectorXd zipped = stan::math::zip_index(delta, loc, time);

  EXPECT_EQ(zipped.size(), loc.size());
  EXPECT_DOUBLE_EQ(zipped[0], delta(0, 0));
  EXPECT_DOUBLE_EQ(zipped[1], delta(0, 0));
  EXPECT_DOUBLE_EQ(zipped[2], delta(4, 0));
  EXPECT_DOUBLE_EQ(zipped[3], delta(1, 1));
  EXPECT_DOUBLE_EQ(zipped[4], delta(2, 2));
  EXPECT_DOUBLE_EQ(zipped[5], delta(2, 2));
  EXPECT_DOUBLE_EQ(zipped[6], delta(3, 3));
  EXPECT_DOUBLE_EQ(zipped[7], delta(3, 0));
  EXPECT_DOUBLE_EQ(zipped[8], delta(4, 0));
}

TEST(MathFunctions, zip_index_errors) {
  Eigen::MatrixXd delta(5, 4);  // 5 rows, 4 columns
  delta.setZero();

  std::vector<int> rows = {1, 2, 3};
  std::vector<int> cols = {1, 2, 3};

  // valid input does not throw
  EXPECT_NO_THROW(stan::math::zip_index(delta, rows, cols));

  // index vectors of different sizes
  std::vector<int> cols_short = {1, 2};
  EXPECT_THROW(stan::math::zip_index(delta, rows, cols_short),
               std::invalid_argument);

  // row index too large (matrix has 5 rows)
  std::vector<int> rows_big = {1, 6, 3};
  EXPECT_THROW(stan::math::zip_index(delta, rows_big, cols), std::out_of_range);

  // row index too small (indices are 1-based)
  std::vector<int> rows_zero = {0, 2, 3};
  EXPECT_THROW(stan::math::zip_index(delta, rows_zero, cols),
               std::out_of_range);

  // column index too large (matrix has 4 columns)
  std::vector<int> cols_big = {1, 5, 3};
  EXPECT_THROW(stan::math::zip_index(delta, rows, cols_big), std::out_of_range);

  // column index too small
  std::vector<int> cols_neg = {1, -1, 3};
  EXPECT_THROW(stan::math::zip_index(delta, rows, cols_neg), std::out_of_range);
}
