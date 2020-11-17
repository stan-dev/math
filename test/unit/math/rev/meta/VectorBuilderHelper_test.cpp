#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraitsRevScal, VectorBuilderHelper_false_true) {
  using stan::VectorBuilderHelper;
  using stan::math::size;
  using stan::math::var;

  var a_var(1);

  VectorBuilderHelper<double, false, true> dvv1(size(a_var));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);
}

TEST(MetaTraitsRevArr, VectorBuilderHelper_false_true) {
  using stan::VectorBuilderHelper;
  using stan::math::size;
  using stan::math::var;
  using std::vector;

  std::vector<var> a_std_vector(3);

  VectorBuilderHelper<double, false, true> dvv2(size(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraitsRevArr, VectorBuilderHelper_true_true) {
  using stan::VectorBuilderHelper;
  using stan::math::size;
  using stan::math::var;
  using std::vector;

  var a_var(1);
  std::vector<var> a_std_vector(3);
  VectorBuilderHelper<double, true, true> dvv1(size(a_var));
  dvv1[0] = 0.0;
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  std::vector<double> data1;
  EXPECT_NO_THROW(data1 = dvv1.data());
  EXPECT_EQ(size(a_var), data1.size());

  VectorBuilderHelper<double, true, true> dvv2(size(a_std_vector));
  dvv2[0] = 0.0;
  dvv2[1] = 1.0;
  dvv2[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(1.0, dvv2[1]);
  EXPECT_FLOAT_EQ(2.0, dvv2[2]);
  std::vector<double> data2;
  EXPECT_NO_THROW(data2 = dvv2.data());
  EXPECT_EQ(size(a_std_vector), data2.size());
}

TEST(MetaTraitsRevMat, VectorBuilderHelper_false_true) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilderHelper;
  using stan::math::size;
  using stan::math::var;

  Matrix<var, Dynamic, 1> a_vector(4);
  Matrix<var, 1, Dynamic> a_row_vector(5);

  VectorBuilderHelper<double, false, true> dvv3(size(a_vector));
  EXPECT_THROW(dvv3[0], std::logic_error);
  EXPECT_THROW(dvv3.data(), std::logic_error);

  VectorBuilderHelper<double, false, true> dvv4(size(a_row_vector));
  EXPECT_THROW(dvv4[0], std::logic_error);
  EXPECT_THROW(dvv3.data(), std::logic_error);
}

TEST(MetaTraitsRevMat, VectorBuilderHelper_true_true) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilderHelper;
  using stan::math::size;
  using stan::math::var;

  Matrix<var, Dynamic, 1> a_vector(4);
  Matrix<var, 1, Dynamic> a_row_vector(5);

  VectorBuilderHelper<double, true, true> dvv3(size(a_vector));
  dvv3[0] = 0.0;
  dvv3[1] = 1.0;
  dvv3[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv3[0]);
  EXPECT_FLOAT_EQ(1.0, dvv3[1]);
  EXPECT_FLOAT_EQ(2.0, dvv3[2]);
  std::vector<double> data3;
  EXPECT_NO_THROW(data3 = dvv3.data());
  EXPECT_EQ(size(a_vector), data3.size());

  VectorBuilderHelper<double, true, true> dvv4(size(a_row_vector));
  dvv4[0] = 0.0;
  dvv4[1] = 1.0;
  dvv4[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv4[0]);
  EXPECT_FLOAT_EQ(1.0, dvv4[1]);
  EXPECT_FLOAT_EQ(2.0, dvv4[2]);
  std::vector<double> data4;
  EXPECT_NO_THROW(data4 = dvv4.data());
  EXPECT_EQ(size(a_row_vector), data4.size());
}
