#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, value_of_rec) {
  using stan::math::value_of_rec;
  double x = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of_rec(x));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(5));
}

TEST(MathMatrixPrimArr, value_of_rec) {
  using stan::math::value_of_rec;
  using std::vector;

  vector<double> a;
  for (size_t i = 0; i < 10; ++i)
    a.push_back(i + 1);

  vector<double> b;
  for (size_t i = 10; i < 15; ++i)
    b.push_back(i + 1);

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  EXPECT_STD_VECTOR_FLOAT_EQ(a, d_a);
  EXPECT_STD_VECTOR_FLOAT_EQ(b, d_b);
}

TEST(MathMatrixPrimMat, value_of_rec) {
  using stan::math::value_of_rec;

  Eigen::Matrix<double, 2, 5> a;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 5; j++)
      a(i, j) = i * 5 + j;
  Eigen::Matrix<double, 5, 1> b;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 1; j++)
      b(i, j) = 10 + i * 5 + j;

  Eigen::MatrixXd d_a = value_of_rec(a);
  Eigen::VectorXd d_b = value_of_rec(b);

  EXPECT_MATRIX_FLOAT_EQ(a, d_a);
  EXPECT_MATRIX_FLOAT_EQ(b, d_b);
}

TEST(MathMatrixPrimMat, value_of_rec_expression) {
  using stan::math::value_of_rec;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(5, 4);
  Eigen::MatrixXd res_a = value_of_rec(2 * a);
  Eigen::MatrixXd correct_a = 2 * a;
  EXPECT_MATRIX_EQ(res_a, correct_a);

  Eigen::VectorXi b = Eigen::VectorXi::Random(7);
  Eigen::VectorXd res_b = value_of_rec(2 * b);
  Eigen::VectorXd correct_b = (2 * b).cast<double>();
  EXPECT_MATRIX_EQ(correct_b, res_b);

  Eigen::ArrayXXd c = a.array();
  Eigen::ArrayXXd res_c = value_of_rec(2 * c);
  Eigen::ArrayXXd correct_c = 2 * c;
  for (int i = 0; i < res_c.size(); i++)
    EXPECT_EQ(correct_c(i), res_c(i));
}

TEST(MathFunctions, value_of_rec_return_type_short_circuit_std_vector) {
  std::vector<double> a(5);
  const std::vector<double> b(5);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(a)),
                            std::vector<double>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(b)),
                            const std::vector<double>&>::value));
}

TEST(MathFunctions, value_of_rec_return_type_short_circuit_vector_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> a(5);
  const Eigen::Matrix<double, Eigen::Dynamic, 1> b(5);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(a)),
                            Eigen::Matrix<double, Eigen::Dynamic, 1>&>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::value_of_rec(b)),
                    const Eigen::Matrix<double, Eigen::Dynamic, 1>&>::value));
}

TEST(MathFunctions, value_of_rec_return_type_short_circuit_row_vector_xd) {
  Eigen::Matrix<double, 1, Eigen::Dynamic> a(5);
  const Eigen::Matrix<double, 1, Eigen::Dynamic> b(5);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(a)),
                            Eigen::Matrix<double, 1, Eigen::Dynamic>&>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::value_of_rec(b)),
                    const Eigen::Matrix<double, 1, Eigen::Dynamic>&>::value));
}

TEST(MathFunctions, value_of_rec_return_type_short_circuit_matrix_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(5, 4);
  const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(5, 4);
  EXPECT_TRUE((std::is_same<
               decltype(stan::math::value_of_rec(a)),
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(b)),
                            const Eigen::Matrix<double, Eigen::Dynamic,
                                                Eigen::Dynamic>&>::value));
}

TEST(MathFunctions, value_of_rec_return_type_short_circuit_expression) {
  const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(5, 4);

  const auto& expr = 3 * a;

  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(expr)),
                            decltype(expr)>::value));
}

TEST(MathFunctions,
     value_of_rec_return_type_short_circuit_static_sized_matrix) {
  Eigen::Matrix<double, 5, 4> a;
  const Eigen::Matrix<double, 5, 4> b;
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(a)),
                            Eigen::Matrix<double, 5, 4>&>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of_rec(b)),
                            const Eigen::Matrix<double, 5, 4>&>::value));
}
