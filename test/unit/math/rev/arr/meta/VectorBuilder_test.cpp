#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorBuilder_false_true) {
  using std::vector;
  using stan::VectorBuilder;
  using stan::math::var;
  using stan::length;

  var a_var(1);
  std::vector<var> a_std_vector(3);

  VectorBuilder<false,double,std::vector<var> > dvv1(length(a_var));
  EXPECT_THROW(dvv1[0], std::logic_error);

  VectorBuilder<false,double,std::vector<var> > dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
}

TEST(MetaTraits, VectorBuilder_true_true) {
  using std::vector;
  using stan::VectorBuilder;
  using stan::math::var;
  using stan::length;

  var a_var(1);
  std::vector<var> a_std_vector(3);

  VectorBuilder<true,double,std::vector<var> > dvv1(length(a_var));
  dvv1[0] = 0.0;
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  
  VectorBuilder<true,double,std::vector<var> > dvv2(length(a_std_vector));
  dvv2[0] = 0.0;
  dvv2[1] = 1.0;
  dvv2[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(1.0, dvv2[1]);
  EXPECT_FLOAT_EQ(2.0, dvv2[2]);  
}

