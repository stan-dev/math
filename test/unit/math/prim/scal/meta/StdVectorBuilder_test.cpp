#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, StdVectorBuilder_false_false) {
  using stan::StdVectorBuilder;
  using stan::length;

  double a_double(1);

  StdVectorBuilder<false, double, double> dvv1(length(a_double));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);
}

TEST(MetaTraits, StdVectorBuilder_true_false) {
  using stan::StdVectorBuilder;
  using stan::length;

  double a_double(1);

  StdVectorBuilder<true, double, double> dvv1(length(a_double));
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  EXPECT_FLOAT_EQ(0.0, dvv1[1]);
  EXPECT_FLOAT_EQ(0.0, dvv1[100]);
  double data1 = 0;
  EXPECT_NO_THROW(data1 = dvv1.data());
  EXPECT_FLOAT_EQ(0.0, data1);
}

TEST(MetaTraits, StdVectorBuilder_type_check) {
  using stan::StdVectorBuilder;
  using stan::contains_std_vector;

  bool r
      = contains_std_vector<StdVectorBuilder<true, double, int>::type>::value;
  EXPECT_FALSE(r);
  r = contains_std_vector<StdVectorBuilder<true, double, double>::type>::value;
  EXPECT_FALSE(r);
}
