#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/sort_test_util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <vector>

template <typename T>
void test_sort_asc() {
  using stan::math::index_type_t;
  using stan::math::sort_asc;

  T c(1);
  c[0] = 1.7;
  T d = sort_asc(c);
  EXPECT_EQ(c.size(), d.size());
  EXPECT_EQ(c[0], d[0]);

  T e(2);
  e[0] = 5.9;
  e[1] = -1.2;
  T f = sort_asc(e);
  EXPECT_EQ(e.size(), f.size());
  EXPECT_EQ(e[0], f[1]);
  EXPECT_EQ(e[1], f[0]);

  T g(3);
  g[0] = 5.9;
  g[1] = -1.2;
  g[2] = 192.13456;
  T h = sort_asc(g);
  EXPECT_EQ(g.size(), h.size());
  EXPECT_EQ(g[0], h[1]);
  EXPECT_EQ(g[1], h[0]);
  EXPECT_EQ(g[2], h[2]);

  T z;
  EXPECT_NO_THROW(sort_asc(z));
  EXPECT_EQ(index_type_t<T>(0), z.size());
}

TEST(MathMatrixPrimMat, sortAscStdVecNan) {
  test_sort_asc_throws<std::vector<double> >();
}

TEST(MathMatrixPrimMat, sortDescStdVecNan) {
  test_sort_desc_throws<std::vector<double> >();
}

TEST(MathMatrixPrimMat, sortAscEigenVecNan) {
  test_sort_asc_throws<Eigen::VectorXd>();
}
TEST(MathMatrixPrimMat, sortAscEigenRowVecNan) {
  test_sort_asc_throws<Eigen::RowVectorXd>();
}

TEST(MathMatrixPrimMat, sortDescEigenVecNan) {
  test_sort_desc_throws<Eigen::VectorXd>();
}
TEST(MathMatrixPrimMat, sortDescEigenRowVecNan) {
  test_sort_desc_throws<Eigen::RowVectorXd>();
}

TEST(MathMatrixPrimMat, sort_asc) {
  using stan::math::sort_asc;

  EXPECT_EQ(0U, sort_asc(std::vector<int>(0)).size());

  test_sort_asc<std::vector<double> >();
  test_sort_asc<Eigen::Matrix<double, Eigen::Dynamic, 1> >();
  test_sort_asc<Eigen::Matrix<double, 1, Eigen::Dynamic> >();
}

template <typename T>
void test_sort_desc() {
  using stan::math::index_type_t;
  using stan::math::sort_desc;

  T c(1);
  c[0] = -1.7;
  T d = sort_desc(c);
  EXPECT_EQ(c.size(), d.size());
  EXPECT_EQ(c[0], d[0]);

  T e(2);
  e[0] = -5.9;
  e[1] = 1.2;
  T f = sort_desc(e);
  EXPECT_EQ(e.size(), f.size());
  EXPECT_EQ(e[0], f[1]);
  EXPECT_EQ(e[1], f[0]);

  T g(3);
  g[0] = -5.9;
  g[1] = 1.2;
  g[2] = -192.13456;
  T h = sort_desc(g);
  EXPECT_EQ(g.size(), h.size());
  EXPECT_EQ(g[0], h[1]);
  EXPECT_EQ(g[1], h[0]);
  EXPECT_EQ(g[2], h[2]);

  T z;
  EXPECT_NO_THROW(sort_desc(z));
  EXPECT_EQ(index_type_t<T>(0), z.size());
}

TEST(MathMatrixPrimMat, sort_desc) {
  using stan::math::sort_desc;

  EXPECT_EQ(0U, sort_desc(std::vector<int>(0)).size());

  test_sort_desc<std::vector<double> >();
  test_sort_desc<Eigen::Matrix<double, Eigen::Dynamic, 1> >();
  test_sort_desc<Eigen::Matrix<double, 1, Eigen::Dynamic> >();
}
