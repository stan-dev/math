#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, VectorBuilder_false_false_scalar) {
  using stan::VectorBuilder;
  using stan::math::size;

  double a_double(1);

  VectorBuilder<false, double, double> dvv1(stan::math::size(a_double));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);
}

TEST(MathMetaPrim, VectorBuilder_true_false_scalar) {
  using stan::VectorBuilder;
  using stan::math::size;

  double a_double(1);

  VectorBuilder<true, double, double> dvv1(stan::math::size(a_double));
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  EXPECT_FLOAT_EQ(0.0, dvv1[1]);
  EXPECT_FLOAT_EQ(0.0, dvv1[100]);
  double data1 = 0;
  EXPECT_NO_THROW(data1 = dvv1.data());
  EXPECT_FLOAT_EQ(0.0, data1);
}

TEST(MathMetaPrim, VectorBuilder_false_false_vector) {
  using stan::VectorBuilder;
  using stan::math::size;
  using std::vector;

  std::vector<double> a_std_vector(3);

  VectorBuilder<false, double, double> dvv2(stan::math::size(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MathMetaPrim, VectorBuilder_true_false_vector) {
  using stan::VectorBuilder;
  using stan::math::size;
  using std::vector;

  std::vector<double> a_std_vector(3);

  VectorBuilder<true, double, double> dvv2(stan::math::size(a_std_vector));
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(0.0, dvv2[1]);
  EXPECT_FLOAT_EQ(0.0, dvv2[2]);
  double data2 = 0;
  EXPECT_NO_THROW(data2 = dvv2.data());
  EXPECT_FLOAT_EQ(0.0, data2);
}

TEST(MathMetaPrim, VectorBuilder_false_false_matrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilder;
  using stan::math::size;

  Matrix<double, Dynamic, 1> a_vector(4);
  Matrix<double, 1, Dynamic> a_row_vector(5);

  VectorBuilder<false, double, double> dvv3(stan::math::size(a_vector));
  EXPECT_THROW(dvv3[0], std::logic_error);
  EXPECT_THROW(dvv3.data(), std::logic_error);

  VectorBuilder<false, double, double> dvv4(stan::math::size(a_row_vector));
  EXPECT_THROW(dvv4[0], std::logic_error);
  EXPECT_THROW(dvv4.data(), std::logic_error);
}

TEST(MathMetaPrim, VectorBuilder_true_false_matrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilder;
  using stan::math::size;

  Matrix<double, Dynamic, 1> a_vector(4);
  Matrix<double, 1, Dynamic> a_row_vector(5);

  VectorBuilder<true, double, double> dvv3(stan::math::size(a_vector));
  EXPECT_FLOAT_EQ(0.0, dvv3[0]);
  EXPECT_FLOAT_EQ(0.0, dvv3[1]);
  EXPECT_FLOAT_EQ(0.0, dvv3[2]);
  double data3 = 0.0;
  EXPECT_NO_THROW(data3 = dvv3.data());
  EXPECT_FLOAT_EQ(0.0, data3);

  VectorBuilder<true, double, double> dvv4(stan::math::size(a_row_vector));
  EXPECT_FLOAT_EQ(0.0, dvv4[0]);
  EXPECT_FLOAT_EQ(0.0, dvv4[1]);
  EXPECT_FLOAT_EQ(0.0, dvv4[2]);
  double data4 = 0.0;
  EXPECT_NO_THROW(data4 = dvv4.data());
  EXPECT_FLOAT_EQ(0.0, data4);
}
