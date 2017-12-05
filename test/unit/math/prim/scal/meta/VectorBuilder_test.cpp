#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorBuilder_false_false) {
  using stan::VectorBuilder;
  using stan::length;

  double a_double(1);

  VectorBuilder<false, double, double> dvv1(length(a_double));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);
}

TEST(MetaTraits, VectorBuilder_true_false) {
  using stan::VectorBuilder;
  using stan::length;

  double a_double(1);

  VectorBuilder<true, double, double> dvv1(length(a_double));
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  EXPECT_FLOAT_EQ(0.0, dvv1[1]);
  EXPECT_FLOAT_EQ(0.0, dvv1[100]);
  double data1 = 0;
  EXPECT_NO_THROW(data1 = dvv1.data());
  EXPECT_FLOAT_EQ(0.0, data1);
}
