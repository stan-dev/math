#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(matrixTest, seq_view_double_stdvector) {
  using stan::math::seq_view;
  std::vector<double> x(2);
  x[0] = 1.0;
  x[1] = 2.0;
  seq_view<double, std::vector<double> > view_x(x);
  EXPECT_EQ(2, view_x.size());
  EXPECT_FLOAT_EQ(1.0, view_x[0]);
  EXPECT_FLOAT_EQ(2.0, view_x[1]);
}

TEST(matrixTest, seq_view_double_vector) {
  using stan::math::seq_view;
  using stan::math::vector_d;
  vector_d y(2);
  y[0] = 1.0;
  y[1] = 2.0;
  seq_view<double, vector_d> view_y(y);
  EXPECT_EQ(2, view_y.size());
  EXPECT_FLOAT_EQ(1.0, view_y[0]);
  EXPECT_FLOAT_EQ(2.0, view_y[1]);
}

TEST(matrixTest, seq_view_double_row_vector) {
  using stan::math::row_vector_d;
  using stan::math::seq_view;
  row_vector_d y(2);
  y[0] = 1.0;
  y[1] = 2.0;
  seq_view<double, row_vector_d> view_y(y);
  EXPECT_EQ(2, view_y.size());
  EXPECT_FLOAT_EQ(1.0, view_y[0]);
  EXPECT_FLOAT_EQ(2.0, view_y[1]);
}

TEST(matrixTest, seq_view_double_double) {
  using stan::math::seq_view;
  double x = 2.0;
  seq_view<double, double> view_x(x);
  EXPECT_EQ(1, view_x.size());
  EXPECT_FLOAT_EQ(2.0, view_x[0]);
  EXPECT_FLOAT_EQ(2.0, view_x[1]);
}

TEST(matrixTest, seq_view_double_matrix) {
  using stan::math::seq_view;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;
  seq_view<double, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
      view_m(m);
  EXPECT_EQ(6, view_m.size());
  for (int i = 0; i < 6; ++i)
    EXPECT_FLOAT_EQ(i + 1, view_m[i]);
}

TEST(matrixTest, seq_view_vec_vec_double) {
  using stan::math::seq_view;
  std::vector<std::vector<double> > x(2);
  for (size_t m = 0; m < 2; ++m)
    x[m] = std::vector<double>(3);
  int pos = 1;
  for (size_t m = 0; m < 2; ++m)
    for (size_t n = 0; n < 3; ++n)
      x[m][n] = pos++;
  seq_view<double, std::vector<std::vector<double> > > view_x(x);
  EXPECT_EQ(6, view_x.size());
  for (int i = 0; i < 6; ++i)
    EXPECT_FLOAT_EQ(i + 1, view_x[i]);
}

TEST(matrixTest, seq_view_double_int) {
  using stan::math::seq_view;
  std::vector<int> x(3);
  x[0] = 0;
  x[1] = 1;
  x[2] = 4;
  seq_view<double, std::vector<int> > view_x(x);
}
