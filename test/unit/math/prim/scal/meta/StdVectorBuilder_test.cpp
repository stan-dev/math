#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraitsPrimScal, StdVectorBuilder_false_false) {
  using stan::StdVectorBuilder;
  using stan::length;

  double a_double(1);

  StdVectorBuilder<false, double, double> dvv1(length(a_double));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);
}

TEST(MetaTraitsPrimScal, StdVectorBuilder_true_false) {
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

TEST(MetaTraitsPrimScal, StdVectorBuilder_type_check) {
  using stan::StdVectorBuilder;
  using stan::contains_std_vector;

  bool r
      = contains_std_vector<StdVectorBuilder<true, double, int>::type>::value;
  EXPECT_FALSE(r);
  r = contains_std_vector<StdVectorBuilder<true, double, double>::type>::value;
  EXPECT_FALSE(r);
}

TEST(MetaTraitsPrimArr, StdVectorBuilder_false_false) {
  using stan::StdVectorBuilder;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  StdVectorBuilder<false, double, double> dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraitsPrimArr, StdVectorBuilder_true_false) {
  using stan::StdVectorBuilder;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  StdVectorBuilder<true, double, double> dvv2(length(a_std_vector));
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(0.0, dvv2[1]);
  EXPECT_FLOAT_EQ(0.0, dvv2[2]);
  double data2 = 0;
  EXPECT_NO_THROW(data2 = dvv2.data());
  EXPECT_FLOAT_EQ(0.0, data2);
}

TEST(MetaTraitsPrimArr, StdVectorBuilder_type_check) {
  using stan::StdVectorBuilder;
  using stan::contains_std_vector;

  bool r = contains_std_vector<
      StdVectorBuilder<true, double, std::vector<int>>::type>::value;
  EXPECT_TRUE(r);
  r = contains_std_vector<
      StdVectorBuilder<true, double, std::vector<double>>::type>::value;
  EXPECT_TRUE(r);
}
