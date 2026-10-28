#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST_F(AgradRev, zip_index_values_and_gradients) {
  Eigen::MatrixXd values(3, 2);
  values << 11, 12, 21, 22, 31, 32;
  stan::math::var_value<Eigen::MatrixXd> x = values;
  const std::vector<int> rows{3, 1, 3, 2, 1};
  const std::vector<int> cols{2, 1, 2, 1, 2};
  const auto y = stan::math::zip_index(x, rows, cols);

  EXPECT_EQ(y.rows(), 5);
  EXPECT_EQ(y.cols(), 1);
  for (int i = 0; i < y.size(); ++i) {
    EXPECT_DOUBLE_EQ(y.val()(i), values(rows[i] - 1, cols[i] - 1));
  }

  Eigen::VectorXd weights(5);
  weights << 2, -3, 5, 7, -11;
  stan::math::var objective
      = stan::math::dot_product(y, weights) + stan::math::sum(x);
  objective.grad();
  Eigen::MatrixXd expected(3, 2);
  expected << -2, -10, 8, 1, 1, 8;
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      EXPECT_DOUBLE_EQ(x.adj()(i, j), expected(i, j));
    }
  }
}

TEST_F(AgradRev, zip_index_lifetime) {
  stan::math::var_value<Eigen::MatrixXd> x = Eigen::MatrixXd::Ones(2, 3);
  auto y = [&x]() {
    std::vector<int> rows{1, 2, 1};
    std::vector<int> cols{3, 1, 3};
    auto result = stan::math::zip_index(x, std::move(rows), std::move(cols));
    rows.assign(3, 2);
    cols.assign(3, 2);
    return result;
  }();
  auto z = stan::math::zip_index(x, std::vector<int>{2}, std::vector<int>{2});
  stan::math::var objective = stan::math::sum(y) + stan::math::sum(z);
  objective.grad();
  EXPECT_DOUBLE_EQ(x.adj()(0, 0), 0);
  EXPECT_DOUBLE_EQ(x.adj()(0, 1), 0);
  EXPECT_DOUBLE_EQ(x.adj()(0, 2), 2);
  EXPECT_DOUBLE_EQ(x.adj()(1, 0), 1);
  EXPECT_DOUBLE_EQ(x.adj()(1, 1), 1);
  EXPECT_DOUBLE_EQ(x.adj()(1, 2), 0);
}

TEST_F(AgradRev, zip_index_matrix_view) {
  Eigen::MatrixXd values(3, 3);
  values << 11, 12, 13, 21, 22, 23, 31, 32, 33;
  stan::math::var_value<Eigen::MatrixXd> x = values;
  auto y = stan::math::zip_index(x.block(1, 1, 2, 2), std::vector<int>{2, 1, 2},
                                 std::vector<int>{1, 2, 1});
  EXPECT_DOUBLE_EQ(y.val()(0), 32);
  EXPECT_DOUBLE_EQ(y.val()(1), 23);
  EXPECT_DOUBLE_EQ(y.val()(2), 32);
  stan::math::sum(y).grad();
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      double expected = i == 2 && j == 1 ? 2 : (i == 1 && j == 2 ? 1 : 0);
      EXPECT_DOUBLE_EQ(x.adj()(i, j), expected);
    }
  }
}

TEST_F(AgradRev, zip_index_row_major) {
  using matrix_t
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  matrix_t values(2, 3);
  values << 11, 12, 13, 21, 22, 23;
  stan::math::var_value<matrix_t> x = values;
  auto y = stan::math::zip_index(x, std::vector<int>{2, 1, 2},
                                 std::vector<int>{1, 3, 1});
  EXPECT_DOUBLE_EQ(y.val()(0), 21);
  EXPECT_DOUBLE_EQ(y.val()(1), 13);
  EXPECT_DOUBLE_EQ(y.val()(2), 21);
  stan::math::sum(y).grad();
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      double expected = i == 1 && j == 0 ? 2 : (i == 0 && j == 2 ? 1 : 0);
      EXPECT_DOUBLE_EQ(x.adj()(i, j), expected);
    }
  }
}

TEST_F(AgradRev, zip_index_empty) {
  const std::vector<int> indices;
  for (int rows : {0, 2}) {
    for (int cols : {0, 3}) {
      stan::math::var_value<Eigen::MatrixXd> x
          = Eigen::MatrixXd::Ones(rows, cols);
      auto y = stan::math::zip_index(x, indices, indices);
      EXPECT_EQ(y.rows(), 0);
      EXPECT_EQ(y.cols(), 1);
      stan::math::sum(y).grad();
      EXPECT_TRUE(x.adj().isZero());
    }
  }
}

TEST_F(AgradRev, zip_index_errors) {
  stan::math::var_value<Eigen::MatrixXd> x = Eigen::MatrixXd::Zero(3, 2);
  const std::vector<int> valid{1, 2};
  const std::vector<int> empty;
  EXPECT_THROW(stan::math::zip_index(x, valid, empty), std::invalid_argument);
  EXPECT_THROW(stan::math::zip_index(x, empty, valid), std::invalid_argument);
  for (int row : {-1, 0, 4}) {
    EXPECT_THROW(stan::math::zip_index(x, std::vector<int>{1, row}, valid),
                 std::out_of_range);
  }
  for (int col : {-1, 0, 3}) {
    EXPECT_THROW(stan::math::zip_index(x, valid, std::vector<int>{1, col}),
                 std::out_of_range);
  }
  for (int rows : {0, 2}) {
    for (int cols : {0, 3}) {
      if (rows != 0 && cols != 0) {
        continue;
      }
      stan::math::var_value<Eigen::MatrixXd> empty_x
          = Eigen::MatrixXd::Zero(rows, cols);
      EXPECT_THROW(stan::math::zip_index(empty_x, valid, valid),
                   std::out_of_range);
    }
  }
}
