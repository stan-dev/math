#include <stan/math/prim/mat.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MathFunctions, adapt_double) {
  double x = 1.0;

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(1, a.size());
  EXPECT_FLOAT_EQ(x, a(0));

  a(0) = 5.0;
  EXPECT_FLOAT_EQ(5.0, a(0));
}

TEST(MathFunctions, adapt_int) {
  int x = 1;

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(0, a.size());
}

TEST(MathFunctions, adapt_std_vector_double) {
  std::vector<double> x = {{1.0, 2.0}};

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(x.size(), a.size());
  for (size_t i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(x[i], a(i));
}

TEST(MathFunctions, adapt_std_vector_int) {
  std::vector<int> x = {{1, 2}};

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(0, a.size());
}

TEST(MathFunctions, adapt_std_eigen_vector_double) {
  Eigen::VectorXd x(2);

  x << 1.0, 2.0;

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(x.size(), a.size());
  for (int i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(x(i), a(i));
}

TEST(MathFunctions, adapt_std_eigen_rowvector_double) {
  Eigen::RowVectorXd x(2);

  x << 1.0, 2.0;

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(x.size(), a.size());
  for (int i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(x(i), a(i));
}

TEST(MathFunctions, adapt_std_eigen_matrix_double) {
  Eigen::MatrixXd x(2, 3);

  x << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(x.size(), a.size());
  for (int i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(x(i), a(i));
}

TEST(MathFunctions, adapt_all_types) {
  double xd = 1.0;
  int xi = 1;
  std::vector<double> xdv = {{1.0, 2.0}};
  std::vector<int> xiv = {{1, 2}};
  Eigen::VectorXd xev(2);
  Eigen::RowVectorXd xrev(2);
  Eigen::MatrixXd xem(2, 3);

  xev << 1.0, 2.0;
  xrev << 2.0, 3.0;
  xem << 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;

  auto a = stan::math::make_variable_adapter<double>(xd, xi, xdv, xiv, xd, xev,
                                                     xrev, xem);

  EXPECT_EQ(2 + xdv.size() + xev.size() + xrev.size() + xem.size(), a.size());
  EXPECT_FLOAT_EQ(xd, a(0));
  for (size_t i = 0; i < xdv.size(); ++i)
    EXPECT_FLOAT_EQ(xdv[i], a(1 + i));
  EXPECT_FLOAT_EQ(xd, a(1 + xdv.size()));

  a(xdv.size()) = 5.0;
  EXPECT_FLOAT_EQ(5.0, a(xdv.size()));
  a(1 + xdv.size()) = 4.0;
  EXPECT_FLOAT_EQ(4.0, a(1 + xdv.size()));
  for (size_t i = 0; i < xev.size(); ++i)
    EXPECT_FLOAT_EQ(xev(i), a(2 + xdv.size() + i));
  for (size_t i = 0; i < xrev.size(); ++i)
    EXPECT_FLOAT_EQ(xrev(i), a(2 + xev.size() + xdv.size() + i));
  for (size_t i = 0; i < xem.size(); ++i)
    EXPECT_FLOAT_EQ(xem(i), a(2 + xrev.size() + xev.size() + xdv.size() + i));
}
