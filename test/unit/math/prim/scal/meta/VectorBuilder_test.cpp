#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/type_traits.hpp>

TEST(MetaTraits, VectorBuilderHelper_false_false) {
  using stan::VectorBuilderHelper;
  using stan::length;

  double a_double(1);

  VectorBuilderHelper<double,false,false> dvv1(length(a_double));
  EXPECT_THROW(dvv1[0], std::logic_error);
}

TEST(MetaTraits, VectorBuilderHelper_true_false) {
  using stan::VectorBuilderHelper;
  using stan::length;

  double a_double(1);

  VectorBuilderHelper<double,true,false> dvv1(length(a_double));
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  EXPECT_FLOAT_EQ(0.0, dvv1[1]);
  EXPECT_FLOAT_EQ(0.0, dvv1[100]);
}

TEST(MetaTraits, VectorBuilder_false_false) {
  using stan::VectorBuilder;
  using stan::length;

  double a_double(1);

  VectorBuilder<false,double,double> dvv1(length(a_double));
  EXPECT_THROW(dvv1[0], std::logic_error);
}

TEST(MetaTraits, VectorBuilder_true_false) {
  using stan::VectorBuilder;
  using stan::length;

  double a_double(1);

  VectorBuilder<true,double,double> dvv1(length(a_double));
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  EXPECT_FLOAT_EQ(0.0, dvv1[1]);
  EXPECT_FLOAT_EQ(0.0, dvv1[100]);
}
