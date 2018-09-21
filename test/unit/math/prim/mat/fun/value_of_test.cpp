#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <type_traits>

TEST(MathMatrix, value_of) {
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

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b(i), d_b(i));

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 5; ++j)
      EXPECT_FLOAT_EQ(a(i, j), d_a(i, j));
}

TEST(MathFunctions, value_of_return_type_short_circuit_vector_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> a(5);
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, Eigen::Dynamic, 1>>::value));
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    Eigen::Matrix<double, Eigen::Dynamic, 1>&>::value));
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    const Eigen::Matrix<double, Eigen::Dynamic, 1>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    const Eigen::Matrix<double, Eigen::Dynamic, 1>&>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_row_vector_xd) {
  Eigen::Matrix<double, 1, Eigen::Dynamic> a(5);
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, 1, Eigen::Dynamic>>::value));
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    Eigen::Matrix<double, 1, Eigen::Dynamic>&>::value));
  EXPECT_FALSE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    const Eigen::Matrix<double, 1, Eigen::Dynamic>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::value_of(a)),
                    const Eigen::Matrix<double, 1, Eigen::Dynamic>&>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_matrix_xd) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(5, 4);
  EXPECT_FALSE((std::is_same<
                decltype(stan::math::value_of(a)),
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>::value));
  EXPECT_FALSE(
      (std::is_same<
          decltype(stan::math::value_of(a)),
          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&>::value));
  EXPECT_FALSE(
      (std::is_same<
          decltype(stan::math::value_of(a)),
          const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            const Eigen::Matrix<double, Eigen::Dynamic,
                                                Eigen::Dynamic>&>::value));
}

TEST(MathFunctions, value_of_return_type_short_circuit_static_sized_matrix) {
  Eigen::Matrix<double, 5, 4> a;
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, 5, 4>>::value));
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             Eigen::Matrix<double, 5, 4>&>::value));
  EXPECT_FALSE((std::is_same<decltype(stan::math::value_of(a)),
                             const Eigen::Matrix<double, 5, 4>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(a)),
                            const Eigen::Matrix<double, 5, 4>&>::value));
}
