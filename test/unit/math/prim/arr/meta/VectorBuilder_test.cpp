#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, VectorBuilder_false_false) {
  using stan::VectorBuilder;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  VectorBuilder<false, double, double> dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraits, VectorBuilder_true_false) {
  using stan::VectorBuilder;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  VectorBuilder<true, double, double> dvv2(length(a_std_vector));
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(0.0, dvv2[1]);
  EXPECT_FLOAT_EQ(0.0, dvv2[2]);
  double data2 = 0;
  EXPECT_NO_THROW(data2 = dvv2.data());
  EXPECT_FLOAT_EQ(0.0, data2);
}
