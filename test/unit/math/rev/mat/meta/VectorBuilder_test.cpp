#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorBuilder_false_true) {
  using stan::VectorBuilder;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;
  using stan::length;

  Matrix<var,Dynamic,1> a_vector(4);
  Matrix<var,1,Dynamic> a_row_vector(5);

  VectorBuilder<false,double,Matrix<var,Dynamic,1> > dvv3(length(a_vector));
  EXPECT_THROW(dvv3[0], std::logic_error);
  EXPECT_THROW(dvv3.data(), std::logic_error);
  
  VectorBuilder<false,double,Matrix<var,1,Dynamic> > dvv4(length(a_row_vector));
  EXPECT_THROW(dvv4[0], std::logic_error);
  EXPECT_THROW(dvv4.data(), std::logic_error);
}

TEST(MetaTraits, VectorBuilder_true_true) {
  using stan::VectorBuilder;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;
  using stan::length;

  Matrix<var,Dynamic,1> a_vector(4);
  Matrix<var,1,Dynamic> a_row_vector(5);

  VectorBuilder<true,double,Matrix<var,Dynamic,1> > dvv3(length(a_vector));
  dvv3[0] = 0.0;
  dvv3[1] = 1.0;
  dvv3[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv3[0]);
  EXPECT_FLOAT_EQ(1.0, dvv3[1]);
  EXPECT_FLOAT_EQ(2.0, dvv3[2]);  
  std::vector<double> data3;
  EXPECT_NO_THROW(data3 = dvv3.data());
  EXPECT_EQ(length(a_vector), data3.size());
  
  VectorBuilder<true,double,Matrix<var,1,Dynamic> > dvv4(length(a_row_vector));
  dvv4[0] = 0.0;
  dvv4[1] = 1.0;
  dvv4[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv4[0]);
  EXPECT_FLOAT_EQ(1.0, dvv4[1]);
  EXPECT_FLOAT_EQ(2.0, dvv4[2]);
  std::vector<double> data4;
  EXPECT_NO_THROW(data4 = dvv4.data());
  EXPECT_EQ(length(a_row_vector), data4.size());
}
