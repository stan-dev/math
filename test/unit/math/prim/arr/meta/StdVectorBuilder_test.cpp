#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, StdVectorBuilder_false_false) {
  using stan::StdVectorBuilder;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  StdVectorBuilder<false, double, double> dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraits, StdVectorBuilder_true_false) {
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

TEST(MetaTraits, StdVectorBuilder_type_check) {
  using stan::StdVectorBuilder;
  using stan::contains_std_vector;

  bool r = contains_std_vector<
      StdVectorBuilder<true, double, std::vector<int>>::type>::value;
  EXPECT_TRUE(r);
  r = contains_std_vector<
      StdVectorBuilder<true, double, std::vector<double>>::type>::value;
  EXPECT_TRUE(r);
}
