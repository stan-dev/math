#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

TEST(MathFunctions, value_of) {
  using stan::math::value_of;
  double x = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of(x));
  EXPECT_FLOAT_EQ(5.0, value_of(5));
}

TEST(MathFunctions, value_of_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::value_of(nan)));
}

TEST(MathMatrixPrimArr, value_of) {
  using stan::math::value_of;
  using std::vector;

  vector<double> a;
  for (size_t i = 0; i < 10; ++i)
    a.push_back(i + 1);

  vector<double> b;
  for (size_t i = 10; i < 15; ++i)
    b.push_back(i + 1);

  vector<double> d_a = value_of(a);
  vector<double> d_b = value_of(b);

  EXPECT_STD_VECTOR_FLOAT_EQ(a, d_a);
  EXPECT_STD_VECTOR_FLOAT_EQ(b, d_b);
}

TEST(MathFunctions, value_of_int_return_type_short_circuit) {
  std::vector<int> a(5, 0);
  const std::vector<int> b(5, 0);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            std::vector<int>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(b)),
                            const std::vector<int>&>::value));
}

TEST(MathFunctions, value_of_double_return_type_short_circuit) {
  std::vector<double> a(5, 0);
  const std::vector<double> b(5, 0);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            std::vector<double>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(b)),
                            const std::vector<double>&>::value));
}

TEST(MathMatrixPrimMat, value_of) {
  using stan::math::value_of;

  Eigen::Matrix<double, 2, 5> a;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++)
      a(i, j) = i * 5 + j;
  Eigen::Matrix<double, 5, 1> b;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 1; j++)
      b(i, j) = 10 + i * 5 + j;

  Eigen::MatrixXd d_a = value_of(a);
  Eigen::VectorXd d_b = value_of(b);

  EXPECT_MATRIX_FLOAT_EQ(a, d_a);
  EXPECT_MATRIX_FLOAT_EQ(b, d_b);
}

TEST(MathFunctions, value_of_vector_of_vectors) {
  std::vector<double> a(5, 0);
  const std::vector<double> b(5, 0);
  std::vector<std::vector<double>> va(5, a);
  const std::vector<std::vector<double>> vb(5, b);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(va)),
                            std::vector<std::vector<double>>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(vb)),
                            const std::vector<std::vector<double>>&>::value));

  auto vva = stan::math::value_of(va);
  auto vvb = stan::math::value_of(va);

  for (size_t i = 0; i < va.size(); ++i) {
    for (size_t j = 0; j < va[i].size(); ++j) {
      EXPECT_FLOAT_EQ(vva[i][j], a[j]);
    }
  }

  for (size_t i = 0; i < vb.size(); ++i) {
    for (size_t j = 0; j < vb[i].size(); ++j) {
      EXPECT_FLOAT_EQ(vvb[i][j], b[j]);
    }
  }
}

TEST(MathFunctions, value_of_vector_of_eigen) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(5);
  Eigen::RowVectorXd b = Eigen::RowVectorXd::Random(5);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(5, 5);
  std::vector<Eigen::VectorXd> va(5, a);
  std::vector<Eigen::RowVectorXd> vb(5, b);
  std::vector<Eigen::MatrixXd> vc(5, c);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(va)),
                            std::vector<Eigen::VectorXd>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(vb)),
                            std::vector<Eigen::RowVectorXd>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(vc)),
                            std::vector<Eigen::MatrixXd>&>::value));

  auto vva = stan::math::value_of(va);
  auto vvb = stan::math::value_of(vb);
  auto vvc = stan::math::value_of(vc);

  for (size_t i = 0; i < vva.size(); ++i)
    EXPECT_MATRIX_FLOAT_EQ(vva[i], a);

  for (size_t i = 0; i < vvb.size(); ++i)
    EXPECT_MATRIX_FLOAT_EQ(vvb[i], b);

  for (size_t i = 0; i < vvc.size(); ++i)
    EXPECT_MATRIX_FLOAT_EQ(vvc[i], c);
}

TEST(MathMatrixPrimMat, value_of_expression) {
  using stan::math::value_of;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(5, 4);
  Eigen::MatrixXd res_a = value_of(2 * a);
  Eigen::MatrixXd correct_a = 2 * a;
  EXPECT_MATRIX_EQ(correct_a, res_a);

  Eigen::VectorXi b = Eigen::VectorXi::Random(7);
  Eigen::VectorXi res_b = value_of(2 * b);
  Eigen::VectorXi correct_b = 2 * b;
  EXPECT_MATRIX_EQ(correct_b, res_b);

  Eigen::ArrayXXd c = a.array();
  Eigen::ArrayXXd res_c = value_of(2 * c);
  Eigen::ArrayXXd correct_c = 2 * c;
  for (int i = 0; i < res_c.size(); i++)
    EXPECT_EQ(correct_c(i), res_c(i));
}

TEST(MathFunctions, value_of_return_type_short_circuit_std_vector) {
  std::vector<double> a(5);
  const std::vector<double> b(5);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            std::vector<double>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(b)),
                            const std::vector<double>&>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_vector_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> a(5);
  const Eigen::Matrix<double, Eigen::Dynamic, 1> b(5);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            Eigen::Matrix<double, Eigen::Dynamic, 1>&>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::value_of(b)),
                    const Eigen::Matrix<double, Eigen::Dynamic, 1>&>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_row_vector_xd) {
  Eigen::Matrix<double, 1, Eigen::Dynamic> a(5);
  const Eigen::Matrix<double, 1, Eigen::Dynamic> b(5);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            Eigen::Matrix<double, 1, Eigen::Dynamic>&>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::value_of(b)),
                    const Eigen::Matrix<double, 1, Eigen::Dynamic>&>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_matrix_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(5, 4);
  const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(5, 4);
  EXPECT_TRUE((std::is_same<
               decltype(stan::math::value_of(a)),
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(b)),
                            const Eigen::Matrix<double, Eigen::Dynamic,
                                                Eigen::Dynamic>&>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_expression) {
  const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(5, 4);

  const auto& expr = 3 * a;

  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(expr)),
                            decltype(expr)>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_static_sized_matrix) {
  Eigen::Matrix<double, 5, 4> a;
  const Eigen::Matrix<double, 5, 4> b;
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            Eigen::Matrix<double, 5, 4>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(b)),
                            const Eigen::Matrix<double, 5, 4>&>::value));
}
