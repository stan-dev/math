#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraitsRevArr, VectorBuilder_false_true) {
  using stan::VectorBuilder;
  using stan::math::size;
  using stan::math::var;
  using std::vector;

  var a_var(1);
  std::vector<var> a_std_vector(3);

  VectorBuilder<false, double, std::vector<var> > dvv1(stan::math::size(a_var));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);

  VectorBuilder<false, double, std::vector<var> > dvv2(
      stan::math::size(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraitsRevArr, VectorBuilder_true_true) {
  using stan::VectorBuilder;
  using stan::math::size;
  using stan::math::var;
  using std::vector;

  var a_var(1);
  std::vector<var> a_std_vector(3);

  VectorBuilder<true, double, std::vector<var> > dvv1(stan::math::size(a_var));
  dvv1[0] = 0.0;
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  std::vector<double> data1;
  EXPECT_NO_THROW(data1 = dvv1.data());
  EXPECT_EQ(stan::math::size(a_var), data1.size());

  VectorBuilder<true, double, std::vector<var> > dvv2(
      stan::math::size(a_std_vector));
  dvv2[0] = 0.0;
  dvv2[1] = 1.0;
  dvv2[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(1.0, dvv2[1]);
  EXPECT_FLOAT_EQ(2.0, dvv2[2]);
  std::vector<double> data2;
  EXPECT_NO_THROW(data2 = dvv2.data());
  EXPECT_EQ(stan::math::size(a_std_vector), data2.size());
}

TEST(MetaTraitsRevMat, VectorBuilder_false_true) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilder;
  using stan::math::size;
  using stan::math::var;

  Matrix<var, Dynamic, 1> a_vector(4);
  Matrix<var, 1, Dynamic> a_row_vector(5);

  VectorBuilder<false, double, Matrix<var, Dynamic, 1> > dvv3(
      stan::math::size(a_vector));
  EXPECT_THROW(dvv3[0], std::logic_error);
  EXPECT_THROW(dvv3.data(), std::logic_error);

  VectorBuilder<false, double, Matrix<var, 1, Dynamic> > dvv4(
      stan::math::size(a_row_vector));
  EXPECT_THROW(dvv4[0], std::logic_error);
  EXPECT_THROW(dvv4.data(), std::logic_error);
}

TEST(MetaTraitsRevMat, VectorBuilder_true_true) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilder;
  using stan::math::size;
  using stan::math::var;

  Matrix<var, Dynamic, 1> a_vector(4);
  Matrix<var, 1, Dynamic> a_row_vector(5);

  VectorBuilder<true, double, Matrix<var, Dynamic, 1> > dvv3(
      stan::math::size(a_vector));
  dvv3[0] = 0.0;
  dvv3[1] = 1.0;
  dvv3[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv3[0]);
  EXPECT_FLOAT_EQ(1.0, dvv3[1]);
  EXPECT_FLOAT_EQ(2.0, dvv3[2]);
  std::vector<double> data3;
  EXPECT_NO_THROW(data3 = dvv3.data());
  EXPECT_EQ(stan::math::size(a_vector), data3.size());

  VectorBuilder<true, double, Matrix<var, 1, Dynamic> > dvv4(
      stan::math::size(a_row_vector));
  dvv4[0] = 0.0;
  dvv4[1] = 1.0;
  dvv4[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv4[0]);
  EXPECT_FLOAT_EQ(1.0, dvv4[1]);
  EXPECT_FLOAT_EQ(2.0, dvv4[2]);
  std::vector<double> data4;
  EXPECT_NO_THROW(data4 = dvv4.data());
  EXPECT_EQ(stan::math::size(a_row_vector), data4.size());
}
