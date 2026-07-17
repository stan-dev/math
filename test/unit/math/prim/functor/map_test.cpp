#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <vector>

TEST(MathFunctor, map_square) {
  std::vector<double> x{1.0, 2.0, 3.0};
  auto y = stan::math::map([](double v) { return v * v; }, x);
  ASSERT_EQ(y.size(), 3u);
  EXPECT_FLOAT_EQ(y[0], 1.0);
  EXPECT_FLOAT_EQ(y[1], 4.0);
  EXPECT_FLOAT_EQ(y[2], 9.0);
}

TEST(MathFunctor, map_type_change) {
  std::vector<int> x{1, 2, 3};
  auto y = stan::math::map([](int v) { return v + 0.5; }, x);
  EXPECT_FLOAT_EQ(y[0], 1.5);
  EXPECT_FLOAT_EQ(y[1], 2.5);
  EXPECT_FLOAT_EQ(y[2], 3.5);
}

TEST(MathFunctor, map_empty) {
  std::vector<double> x;
  auto y = stan::math::map([](double v) { return v; }, x);
  EXPECT_EQ(y.size(), 0u);
}

TEST(MathFunctor, map_const_vector) {
  const std::vector<double> x{1.0, 2.0, 3.0};
  auto y = stan::math::map([](double v) { return 2.0 * v; }, x);
  ASSERT_EQ(y.size(), 3u);
  EXPECT_FLOAT_EQ(y[0], 2.0);
  EXPECT_FLOAT_EQ(y[1], 4.0);
  EXPECT_FLOAT_EQ(y[2], 6.0);
}

TEST(MathFunctor, map_variadic_args) {
  std::vector<double> x{1.0, 2.0, 3.0};
  double offset = 10.0;
  double scale = 2.0;
  auto y = stan::math::map(
      [](double v, double mu, double sigma) { return mu + sigma * v; }, x,
      offset, scale);
  ASSERT_EQ(y.size(), 3u);
  EXPECT_FLOAT_EQ(y[0], 12.0);
  EXPECT_FLOAT_EQ(y[1], 14.0);
  EXPECT_FLOAT_EQ(y[2], 16.0);
}

TEST(MathFunctor, map_shared_arg_reused) {
  std::vector<double> x{1.0, 2.0, 3.0};
  auto offset = std::make_unique<double>(10.0);
  auto y = stan::math::map(
      [](double v, const std::unique_ptr<double>& mu) { return v + *mu; }, x,
      std::move(offset));
  ASSERT_EQ(y.size(), 3u);
  EXPECT_FLOAT_EQ(y[0], 11.0);
  EXPECT_FLOAT_EQ(y[1], 12.0);
  EXPECT_FLOAT_EQ(y[2], 13.0);
}

TEST(MathFunctor, mapN_add) {
  std::vector<double> a{1.0, 2.0};
  std::vector<double> b{10.0, 20.0};
  auto y = stan::math::mapN([](double u, double v) { return u + v; }, a, b);
  ASSERT_EQ(y.size(), 2u);
  EXPECT_FLOAT_EQ(y[0], 11.0);
  EXPECT_FLOAT_EQ(y[1], 22.0);
}

TEST(MathFunctor, mapN_mixed_scalar_types) {
  std::vector<int> a{1, 2, 3};
  std::vector<double> b{0.5, 1.5, 2.5};
  auto y = stan::math::mapN([](int u, double v) { return u + v; }, a, b);
  ASSERT_EQ(y.size(), 3u);
  EXPECT_FLOAT_EQ(y[0], 1.5);
  EXPECT_FLOAT_EQ(y[1], 3.5);
  EXPECT_FLOAT_EQ(y[2], 5.5);
}

TEST(MathFunctor, mapN_three_vectors) {
  std::vector<double> a{1.0, 2.0};
  std::vector<double> b{10.0, 20.0};
  std::vector<double> c{100.0, 200.0};
  auto y = stan::math::mapN(
      [](double u, double v, double w) { return u + v + w; }, a, b, c);
  ASSERT_EQ(y.size(), 2u);
  EXPECT_FLOAT_EQ(y[0], 111.0);
  EXPECT_FLOAT_EQ(y[1], 222.0);
}

TEST(MathFunctor, mapN_empty) {
  std::vector<double> a;
  std::vector<double> b;
  auto y = stan::math::mapN([](double u, double v) { return u + v; }, a, b);
  EXPECT_EQ(y.size(), 0u);
}

TEST(MathFunctor, mapN_size_mismatch) {
  std::vector<double> a{1.0, 2.0};
  std::vector<double> b{10.0};
  EXPECT_THROW(stan::math::mapN([](double u, double v) { return u + v; }, a, b),
               std::invalid_argument);
}

TEST(MathFunctor, row_map_standardize) {
  Eigen::MatrixXd x(2, 3);
  x << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  Eigen::RowVectorXd mu(3);
  mu << 1.0, 2.0, 3.0;
  Eigen::RowVectorXd sigma(3);
  sigma << 1.0, 2.0, 3.0;

  auto y = stan::math::row_map(
      [](const Eigen::RowVectorXd& row, const Eigen::RowVectorXd& m,
         const Eigen::RowVectorXd& s) {
        return (row.array() - m.array()) / s.array();
      },
      x, mu, sigma);

  Eigen::MatrixXd expected(2, 3);
  expected << 0.0, 0.0, 0.0, 3.0, 1.5, 1.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, row_map_empty) {
  Eigen::MatrixXd x(0, 3);
  auto y = stan::math::row_map(
      [](const Eigen::RowVectorXd& row) { return row; }, x);
  EXPECT_EQ(y.rows(), 0);
  EXPECT_EQ(y.cols(), 0);
}

TEST(MathFunctor, row_map_inconsistent_row_size) {
  Eigen::MatrixXd x(2, 2);
  x << 1.0, 2.0, 3.0, 4.0;
  int call = 0;
  EXPECT_THROW(stan::math::row_map(
                   [&call](const Eigen::RowVectorXd&) {
                     ++call;
                     if (call == 1) {
                       return Eigen::RowVectorXd::Ones(2);
                     }
                     return Eigen::RowVectorXd::Ones(3);
                   },
                   x),
               std::invalid_argument);
}

TEST(MathFunctor, row_map_shared_arg_reused) {
  Eigen::MatrixXd x(2, 2);
  x << 1.0, 2.0, 3.0, 4.0;
  auto offset = std::make_unique<double>(10.0);
  auto y = stan::math::row_map(
      [](const Eigen::RowVectorXd& row, const std::unique_ptr<double>& mu) {
        return row.array() + *mu;
      },
      x, std::move(offset));

  Eigen::MatrixXd expected(2, 2);
  expected << 11.0, 12.0, 13.0, 14.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, col_map_sum_rows) {
  Eigen::MatrixXd x(2, 2);
  x << 1.0, 2.0, 3.0, 4.0;

  auto y = stan::math::col_map(
      [](const Eigen::VectorXd& col) {
        return Eigen::VectorXd::Constant(1, col.sum());
      },
      x);

  Eigen::MatrixXd expected(1, 2);
  expected << 4.0, 6.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, col_map_empty) {
  Eigen::MatrixXd x(3, 0);
  auto y
      = stan::math::col_map([](const Eigen::VectorXd& col) { return col; }, x);
  EXPECT_EQ(y.rows(), 0);
  EXPECT_EQ(y.cols(), 0);
}

TEST(MathFunctor, col_map_inconsistent_col_size) {
  Eigen::MatrixXd x(2, 2);
  x << 1.0, 2.0, 3.0, 4.0;
  int call = 0;
  EXPECT_THROW(stan::math::col_map(
                   [&call](const Eigen::VectorXd&) {
                     ++call;
                     if (call == 1) {
                       return Eigen::VectorXd::Ones(2);
                     }
                     return Eigen::VectorXd::Ones(3);
                   },
                   x),
               std::invalid_argument);
}

TEST(MathFunctor, col_map_shared_arg_reused) {
  Eigen::MatrixXd x(2, 2);
  x << 1.0, 2.0, 3.0, 4.0;
  auto scale = std::make_unique<double>(2.0);
  auto y = stan::math::col_map(
      [](const Eigen::VectorXd& col, const std::unique_ptr<double>& s) {
        return (*s) * col;
      },
      x, std::move(scale));

  Eigen::MatrixXd expected(2, 2);
  expected << 2.0, 4.0, 6.0, 8.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, row_mapN_add) {
  Eigen::MatrixXd a(2, 2);
  a << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd b(2, 2);
  b << 10.0, 20.0, 30.0, 40.0;

  auto y
      = stan::math::row_mapN([](const Eigen::RowVectorXd& u,
                                const Eigen::RowVectorXd& v) { return u + v; },
                             a, b);

  Eigen::MatrixXd expected(2, 2);
  expected << 11.0, 22.0, 33.0, 44.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, row_mapN_dim_mismatch) {
  Eigen::MatrixXd a(2, 2);
  a << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd b(2, 3);
  b << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  EXPECT_THROW(
      stan::math::row_mapN([](const Eigen::RowVectorXd& u,
                              const Eigen::RowVectorXd& v) { return u + v; },
                           a, b),
      std::invalid_argument);
}

TEST(MathFunctor, row_mapN_empty) {
  Eigen::MatrixXd a(0, 2);
  Eigen::MatrixXd b(0, 2);
  auto y
      = stan::math::row_mapN([](const Eigen::RowVectorXd& u,
                                const Eigen::RowVectorXd& v) { return u + v; },
                             a, b);
  EXPECT_EQ(y.rows(), 0);
  EXPECT_EQ(y.cols(), 0);
}

TEST(MathFunctor, row_mapN_mixed_scalar_types) {
  Eigen::MatrixXi a(2, 2);
  a << 1, 2, 3, 4;
  Eigen::MatrixXd b(2, 2);
  b << 0.5, 1.5, 2.5, 3.5;
  auto y = stan::math::row_mapN(
      [](const auto& u, const auto& v) {
        return u.template cast<double>() + v;
      },
      a, b);

  Eigen::MatrixXd expected(2, 2);
  expected << 1.5, 3.5, 5.5, 7.5;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, row_mapN_inconsistent_row_size) {
  Eigen::MatrixXd a(2, 2);
  a << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd b(2, 2);
  b << 10.0, 20.0, 30.0, 40.0;
  int call = 0;
  EXPECT_THROW(
      stan::math::row_mapN(
          [&call](const Eigen::RowVectorXd&, const Eigen::RowVectorXd&) {
            ++call;
            if (call == 1) {
              return Eigen::RowVectorXd::Ones(2);
            }
            return Eigen::RowVectorXd::Ones(3);
          },
          a, b),
      std::invalid_argument);
}

TEST(MathFunctor, col_mapN_add) {
  Eigen::MatrixXd a(2, 2);
  a << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd b(2, 2);
  b << 10.0, 20.0, 30.0, 40.0;

  auto y = stan::math::col_mapN(
      [](const Eigen::VectorXd& u, const Eigen::VectorXd& v) { return u + v; },
      a, b);

  Eigen::MatrixXd expected(2, 2);
  expected << 11.0, 22.0, 33.0, 44.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, col_mapN_dim_mismatch) {
  Eigen::MatrixXd a(2, 2);
  a << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd b(3, 2);
  b << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  EXPECT_THROW(
      stan::math::col_mapN([](const Eigen::VectorXd& u,
                              const Eigen::VectorXd& v) { return u + v; },
                           a, b),
      std::invalid_argument);
}

TEST(MathFunctor, col_mapN_empty) {
  Eigen::MatrixXd a(2, 0);
  Eigen::MatrixXd b(2, 0);
  auto y = stan::math::col_mapN(
      [](const Eigen::VectorXd& u, const Eigen::VectorXd& v) { return u + v; },
      a, b);
  EXPECT_EQ(y.rows(), 0);
  EXPECT_EQ(y.cols(), 0);
}

TEST(MathFunctor, col_mapN_mixed_scalar_types) {
  Eigen::MatrixXi a(2, 2);
  a << 1, 2, 3, 4;
  Eigen::MatrixXd b(2, 2);
  b << 0.5, 1.5, 2.5, 3.5;
  auto y = stan::math::col_mapN(
      [](const auto& u, const auto& v) {
        return u.template cast<double>() + v;
      },
      a, b);

  Eigen::MatrixXd expected(2, 2);
  expected << 1.5, 3.5, 5.5, 7.5;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, col_mapN_inconsistent_col_size) {
  Eigen::MatrixXd a(2, 2);
  a << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd b(2, 2);
  b << 10.0, 20.0, 30.0, 40.0;
  int call = 0;
  EXPECT_THROW(stan::math::col_mapN(
                   [&call](const Eigen::VectorXd&, const Eigen::VectorXd&) {
                     ++call;
                     if (call == 1) {
                       return Eigen::VectorXd::Ones(2);
                     }
                     return Eigen::VectorXd::Ones(3);
                   },
                   a, b),
               std::invalid_argument);
}

TEST(MathFunctor, row_map_block_expression) {
  Eigen::MatrixXd x(3, 3);
  x << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
  auto y = stan::math::row_map([](const auto& row) { return 2.0 * row; },
                               x.block(0, 0, 2, 3));

  Eigen::MatrixXd expected(2, 3);
  expected << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}

TEST(MathFunctor, col_map_block_expression) {
  Eigen::MatrixXd x(3, 3);
  x << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
  auto y = stan::math::col_map([](const auto& col) { return 2.0 * col; },
                               x.block(0, 0, 3, 2));

  Eigen::MatrixXd expected(3, 2);
  expected << 2.0, 4.0, 8.0, 10.0, 14.0, 16.0;
  EXPECT_MATRIX_FLOAT_EQ(y, expected);
}
