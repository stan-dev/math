#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>


template <typename V, typename T>
bool is_type(const T& x) {
  using stan::math::abs;
  return std::is_same<V, T>::value;
}


TEST(mixFun, absTypes) {
  using stan::math::abs;
  EXPECT_TRUE(is_type<int>(abs(1)));
  EXPECT_TRUE(is_type<double>(abs(-2.3)));
  EXPECT_TRUE(is_type<std::vector<int>>(abs(std::vector<int>{1, 2, -3})));
  EXPECT_TRUE(is_type<std::vector<double>>(abs(std::vector<double>{2.0, 3.9})));
}

TEST(MathFunctions, absInt) {
  using stan::math::abs;

  int x = -3;
  EXPECT_EQ(3, abs(x));
  EXPECT_TRUE(is_type<int>(abs(x)));
}

TEST(MathFunctions, absIntVec) {
  using stan::math::abs;
  using sv_t = std::vector<int>;
  sv_t a;
  EXPECT_EQ(0, abs(a).size());

  sv_t b = { -1, 2 };
  EXPECT_EQ(2, abs(b).size());
  EXPECT_EQ(1, abs(b)[0]);
  EXPECT_EQ(2, abs(b)[1]);

  using svv_t = std::vector<sv_t>;
  svv_t c;
  EXPECT_EQ(0, abs(c).size());

  svv_t d = {{1, -2}, {-101, 15}};
  EXPECT_EQ(2, abs(d).size());
  EXPECT_EQ(2, abs(d)[0].size());
  EXPECT_EQ(1, abs(d)[0][0]);
  EXPECT_EQ(2, abs(d)[0][1]);
  EXPECT_EQ(101, abs(d)[1][0]);
  EXPECT_EQ(15, abs(d)[1][1]);

  using svvv_t = std::vector<svv_t>;
  svvv_t e = {{{-4}}};
  EXPECT_EQ(4, abs(e)[0][0][0]);
}

TEST(MathFunctions, absDouble) {
  double y = 2.0;
  EXPECT_FLOAT_EQ(2.0, abs(y));

  y = 128745.72;
  EXPECT_FLOAT_EQ(128745.72, abs(y));

  y = -y;
  EXPECT_FLOAT_EQ(128745.72, abs(y));

  y = -1.3;
  EXPECT_FLOAT_EQ(1.3, abs(y));
}

TEST(MathFunctions, absDoubleVec) {
  using stan::math::abs;
  using sv_t = std::vector<double>;
  sv_t a;
  EXPECT_EQ(0, abs(a).size());

  sv_t b = { -1.7, 2.9 };
  EXPECT_EQ(2, abs(b).size());
  EXPECT_FLOAT_EQ(1.7, abs(b)[0]);
  EXPECT_FLOAT_EQ(2.9, abs(b)[1]);

  using svv_t = std::vector<sv_t>;
  svv_t c;
  EXPECT_EQ(0, abs(c).size());

  svv_t d = {{1.3, -2}, {-101.27, 15.39}};
  EXPECT_EQ(2, abs(d).size());
  EXPECT_EQ(2, abs(d)[0].size());
  EXPECT_FLOAT_EQ(1.3, abs(d)[0][0]);
  EXPECT_FLOAT_EQ(2, abs(d)[0][1]);
  EXPECT_FLOAT_EQ(101.27, abs(d)[1][0]);
  EXPECT_FLOAT_EQ(15.39, abs(d)[1][1]);
}

TEST(MathFunctions, abs_complex) {
  using c_t = std::complex<double>;
  EXPECT_FLOAT_EQ(5.0, abs(c_t(3.0, -4.0)));
}


TEST(MathFunctions, absComplexVec) {
  using stan::math::abs;
  using c_t = std::complex<double>;
  using sv_t = std::vector<c_t>;
  sv_t a;
  EXPECT_EQ(0, abs(a).size());

  sv_t b = { c_t(3, 4), c_t(-3, 4) };
  EXPECT_EQ(2, abs(b).size());
  EXPECT_FLOAT_EQ(5, abs(b)[0]);
  EXPECT_FLOAT_EQ(5, abs(b)[1]);

  using svv_t = std::vector<sv_t>;
  svv_t c;
  EXPECT_EQ(0, abs(c).size());

  svv_t d = {{1.3, c_t(-std::sqrt(2), -std::sqrt(2))}, {c_t(3, -4), 15.39}};
  EXPECT_EQ(2, abs(d).size());
  EXPECT_EQ(2, abs(d)[0].size());
  EXPECT_FLOAT_EQ(1.3, abs(d)[0][0]);
  EXPECT_FLOAT_EQ(2, abs(d)[0][1]);
  EXPECT_FLOAT_EQ(5, abs(d)[1][0]);
  EXPECT_FLOAT_EQ(15.39, abs(d)[1][1]);
}

TEST(MathFunctions, abs_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::abs(nan)));
  EXPECT_TRUE(std::isnan(stan::math::abs(std::complex<double>(nan, 1.0))));
  EXPECT_TRUE(std::isnan(stan::math::abs(std::complex<double>(nan, nan))));
  EXPECT_TRUE(std::isnan(stan::math::abs(std::complex<double>(1.0, nan))));
}

TEST(MathFunctions, abs_vec) {
  using stan::math::abs;
  Eigen::VectorXd y(0);
  EXPECT_EQ(0, abs(y).size());

  Eigen::VectorXd z(3);
  z << 1.3, -2.9, 0;
  EXPECT_EQ(3, abs(z).size());
  EXPECT_EQ(1.3, abs(z)[0]);
  EXPECT_EQ(2.9, abs(z)[1]);
  EXPECT_EQ(0, abs(z)[2]);
}
TEST(MathFunctions, abs_rowvec) {
  using stan::math::abs;
  Eigen::RowVectorXd y(0);
  EXPECT_EQ(0, abs(y).size());

  Eigen::RowVectorXd z(3);
  z << 1.3, -2.9, 0;
  EXPECT_EQ(3, abs(z).size());
  EXPECT_EQ(1.3, abs(z)[0]);
  EXPECT_EQ(2.9, abs(z)[1]);
  EXPECT_EQ(0, abs(z)[2]);
}
TEST(MathFunctions, abs_mat) {
  using stan::math::abs;
  Eigen::MatrixXd y(0, 0);
  EXPECT_EQ(0, abs(y).size());

  Eigen::MatrixXd z(2, 3);
  z << 1.3, -2.9, 0, 5.8, -112, 0;
  EXPECT_EQ(2, abs(z).rows());
  EXPECT_EQ(3, abs(z).cols());
  EXPECT_FLOAT_EQ(1.3, abs(z)(0, 0));
  EXPECT_FLOAT_EQ(2.9, abs(z)(0, 1));
  EXPECT_FLOAT_EQ(0, abs(z)(0, 2));
  EXPECT_FLOAT_EQ(5.8, abs(z)(1, 0));
  EXPECT_FLOAT_EQ(112, abs(z)(1, 1));
  EXPECT_FLOAT_EQ(0, abs(z)(1, 2));
}

TEST(MathFunctions, abs_complex_vec) {
  using c_t = std::complex<double>;
  using stan::math::abs;
  Eigen::VectorXcd y(0);
  EXPECT_EQ(0, abs(y).size());

  Eigen::VectorXcd z(3);
  z << 1.3, c_t(3, -4), 0;
  EXPECT_EQ(3, abs(z).size());
  EXPECT_EQ(1.3, abs(z)[0]);
  EXPECT_EQ(5, abs(z)[1]);
  EXPECT_EQ(0, abs(z)[2]);
}

TEST(MathFunctions, abs_complex_rowvec) {
  using c_t = std::complex<double>;
  using stan::math::abs;
  Eigen::RowVectorXcd y(0);
  EXPECT_EQ(0, abs(y).size());

  Eigen::RowVectorXcd z(3);
  z << 1.3, c_t(3, -4), 0;
  EXPECT_EQ(3, abs(z).size());
  EXPECT_EQ(1.3, abs(z)[0]);
  EXPECT_EQ(5, abs(z)[1]);
  EXPECT_EQ(0, abs(z)[2]);
}

TEST(MathFunctions, abs_complex_mat) {
  using c_t = std::complex<double>;
  using stan::math::abs;
  Eigen::MatrixXcd y(0, 0);
  EXPECT_EQ(0, abs(y).size());

  Eigen::MatrixXcd z(2, 3);
  z << 1.3, c_t(3, -4), 0, 1, c_t(5, -5), c_t(-1.3, 2.9);
  EXPECT_EQ(2, abs(z).rows());
  EXPECT_EQ(3, abs(z).cols());

  EXPECT_FLOAT_EQ(1.3, abs(z)(0, 0));
  EXPECT_FLOAT_EQ(5, abs(z)(0, 1));
  EXPECT_FLOAT_EQ(0, abs(z)(0, 2));
  EXPECT_FLOAT_EQ(1, abs(z)(1, 0));
  EXPECT_FLOAT_EQ(std::sqrt(50), abs(z)(1, 1));
  EXPECT_FLOAT_EQ(std::sqrt(1.3 * 1.3 + 2.9 * 2.9), abs(z)(1, 2));
}
